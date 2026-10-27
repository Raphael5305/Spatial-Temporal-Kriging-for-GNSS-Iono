function results = STthreat2(ST_STmodels,R, N, IPPVD, samp_int, Rin_numbers,IGPs, Broadcast_GIVD, Broadcast_IGPs,dq)
% Inputs:
% R = [R_min, R_max], N = [N_min, N_target], IPPVD = PF data
% ST_S, ST_T, ST_ST = Data for computing ST_STmodel2
% epoch_t = Epoch time array

if nargin < 10, dq = []; end

R_min = R(1); %#ok<NASGU> % 仅保留 800km for WAAS
R_max = R(2); % 2100km 
N_min = N(1);  %10
N_target = N(2); % 30

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350000/1000; % km
R_earth = 6378.137;    % km

PF = IPPVD; % IPPVD should be detrended PF data:STPF,unit:m


% 规范 IGPs 经度 [-180,180]
IGP_lat = IGPs(:,1);
IGP_lon = mod(IGPs(:,2) + 180, 360) - 180;

% 规范 Broadcast_IGPs 经度
B_lat = Broadcast_IGPs(:,1);
B_lon = mod(Broadcast_IGPs(:,2) + 180, 360) - 180;
Broadcast_IGPs = [B_lat, B_lon];

% 历元集合（按 samp_int 下采样）
all_epochs = unique(IPPVD(:,1), 'stable');
est_epochs = all_epochs(1:samp_int:numel(all_epochs));
num_epochs = numel(est_epochs);
num_igps   = numel(IGP_lat);

save_row=0;
results = zeros((Rin_numbers+1)* num_epochs*num_igps, 16);

% --- IGP → Broadcast_IGPs 映射（精确匹配；你需要最近邻可自行加上） ---
tol = 1e-8;
B_map = nan(num_igps,1);
for ig = 1:num_igps
    mask_lat = abs(B_lat - IGP_lat(ig)) <= tol;
    mask_lon = abs(B_lon - IGP_lon(ig)) <= tol;
    idx = find(mask_lat & mask_lon, 1, 'first');
    if ~isempty(idx)
        B_map(ig) = idx;
    end
end

% ---- 主循环 ----
for ie = 6:num_epochs
    est_time = est_epochs(ie);
    % 该历元的 STmodel2（若不存在或为空 → 走 UnK_threat 分支）
    model_idx = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    has_ST = (~isempty(model_idx) && ~isempty(ST_STmodels(model_idx).model));

    if has_ST
        ST_ST2model = ST_STmodels(model_idx).model;
    else
        continue; % 显式为空用于分支判断
    end

    % 时间窗候选（ST 分支用过去300s；UnK_threat 将改为仅用当前历元）
    time_range = est_time - 300 : est_time;
    time_idx_win   = find(ismember(IPPVD(:,1), time_range));  % 窗口索引
    time_idx_exact = find(IPPVD(:,1) == est_time);            % 当前历元索引

    if isempty(time_idx_win)||isempty(time_idx_exact)
        % ST需要窗口:若窗口空，直接跳过这个历元；若当前无历元也跳过
        continue;
    end

    % 日内 30s 索引 j（列1=00:00:00）
    j = epoch_to_30s_index(est_time);
    jj= epoch_to_30s_index(est_time)-1;%判断哪个时间是对的

    for ig = 1:num_igps
        % Broadcast GIVD 真值（第2列 & 传给 STkriging 用于诊断）
        truth_val = NaN;
        truth_val2= NaN;
        if ~isnan(jj)
            m = B_map(ig);
            if ~isnan(m), truth_val = Broadcast_GIVD(m, j); truth_val2 = Broadcast_GIVD(m, jj);end
        end
        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);

        % 进度提示

        if ~isempty(dq) && mod(ie, 50) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('线程 %d | epoch %d/%d | IGP %d/%d', tid, ie, num_epochs, ig, num_igps));
        end

        % ====================== ST 分支：按 ST 变差值排序 + STkriging ======================
        lat_cand = IPPVD(time_idx_win, 3);
        lon_cand = IPPVD(time_idx_win, 4);

        % 候选样本与 IGP 的大圆距离
        dlat = deg2rad(lat_cand - lat1);
        dlon = deg2rad(lon_cand - lon1);
        a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat_cand)).*sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        Dist_all = R_earth * c;

        selected = (Dist_all <= R_max & Dist_all > 0);
        idsearch_initial = time_idx_win(selected);
        if numel(idsearch_initial) < N_min
            continue;
        end

        % Step1：以 ST 变差函数排序
        s = Dist_all(selected);
        tdiff = abs(IPPVD(idsearch_initial,1) - est_time);
        vvals = ST_ST2model.func(s, tdiff);
        [sorted_variogram, ord] = sort(vvals, 'ascend');%从小到大排列
        idsearch_initial = idsearch_initial(ord);
        variogram_values = sorted_variogram;
        % If N_min <= count <= N_target, use all points
        % Compute max and min variogram values
        max_variogram = max(variogram_values);
        min_variogram = min(variogram_values);
        % Use all IPPs within Rmax to generate GIVD, GIVE
        PF_samples = PF(idsearch_initial, :);%variogram从小到大排列的IPPs
        R_out = 0.5 * (max_variogram+min_variogram);

        R_in_values = linspace(min_variogram-1e-10,max_variogram+1e-10,Rin_numbers+1);
        target_point=[est_time,ig,lat1,lon1,truth_val2];%MT18 as GIVD broadcast

        for r = 1:length(R_in_values)
            R_in = R_in_values(r);
            result_once = STdepriveAndEstimate2(R_min,R_max,N,R_in, R_out, PF_samples, target_point, ST_ST2model, N_min,N_target, wgs84, height, R_earth);
            if isnan(result_once)
                continue;
            end
            save_row = save_row + 1;
            results(save_row, :) = [target_point result_once];
        end
    end
end
% 裁剪
results = results(1:save_row, :);
end

% === 辅助函数：列1=00:00:00 → j = floor(sod/30)+1 ===
function j = epoch_to_30s_index(epoch_seconds)
    if ~isfinite(epoch_seconds), j = NaN; return; end
    sod = mod(epoch_seconds, 86400);
    if sod < 0, sod = sod + 86400; end
    j = floor(sod / 30) + 1; % 0→1, 30→2, ..., 86370→2880
end


 
