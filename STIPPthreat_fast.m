function results = STIPPthreat_fast(ST_STmodels,R, N, IPPVD, samp_int, Rin_numbers,IGPs, Broadcast_GIVD, Broadcast_IGPs,dq)
% 结果严格一致版（仅做等价提速，不改任何选择/顺序/边界判定）
% - 保留容差匹配 Broadcast_IGPs
% - 保留 sort 升序 + 取前 N_target
% - 保留 target_point 使用 truth_val(j)
% - 保留 save_row 逐行写入顺序

if nargin < 10, dq = []; end

R_min    = R(1); % 800
R_max    = R(2); % 2100
N_min    = N(1); % 10
N_target = N(2); % 30

% 常量
wgs84   = wgs84Ellipsoid('kilometer'); %#ok<NASGU>
height  = 350;          % km
R_earth = 6378.137;     % km
DEG2RAD = pi/180;

PF        = IPPVD; % [epoch, sat, lat, lon, ..., y]
pf_epoch  = PF(:,1);
pf_lat    = PF(:,3);
pf_lon    = PF(:,4);

% 规范 IGPs 经度 [-180,180]
IGP_lat = IGPs(:,1);
IGP_lon = mod(IGPs(:,2) + 180, 360) - 180;

% 规范 Broadcast_IGPs 经度（保持与你原逻辑一致）
B_lat = Broadcast_IGPs(:,1);
B_lon = mod(Broadcast_IGPs(:,2) + 180, 360) - 180;

% 历元集合（按 samp_int 下采样）——与原法保持一致
all_epochs = unique(PF(:,1), 'stable');
est_epochs = all_epochs(1:samp_int:numel(all_epochs));
num_epochs = numel(est_epochs);
num_igps   = numel(IGP_lat);

% —— 预分配：保持与你原式同一上界，保证写入顺序与 save_row 一致 —— 
results  = zeros((Rin_numbers*9*5)* num_epochs*num_igps, 16);
save_row = 0;

% --- IGP → Broadcast_IGPs 映射（容差匹配，保持原行为） ---
tol  = 1e-8;
B_map = nan(num_igps,1);
for ig = 1:num_igps
    mask_lat = abs(B_lat - IGP_lat(ig)) <= tol;
    mask_lon = abs(B_lon - IGP_lon(ig)) <= tol;
    idx = find(mask_lat & mask_lon, 1, 'first');
    if ~isempty(idx)
        B_map(ig) = idx;
    end
end

% ---- 主循环（顺序与原来一致：epoch 外层，IGP 内层）----
for ie = 1:num_epochs
    est_time = est_epochs(ie);

    % 与原版一致：仅当该历元存在 STmodel2
    model_idx = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    has_ST = (~isempty(model_idx) && ~isempty(ST_STmodels(model_idx).packed));
    if ~has_ST
        continue;
    end
    ST_ST2model = ST_STmodels(model_idx).packed;

    % 时间窗口（等价于 ismember 的逻辑掩码，不改变顺序）
    mask_win      = (pf_epoch >= est_time-300) & (pf_epoch <= est_time);
    time_idx_win  = find(mask_win);
    time_idx_exact= find(pf_epoch == est_time); %#ok<NASGU> % 与原代码一致地保留
    if isempty(time_idx_win) || isempty(time_idx_exact)
        continue;
    end

    % 日内 30s 索引 j（与原一致）
    j  = epoch_to_30s_index(est_time);
    jj = j - 1; %#ok<NASGU> % 原代码里存在但未用于 target_point，这里同样保留

    for ig = 1:num_igps
        % Broadcast GIVD 真值（保持 j，用于 target_point）
        truth_val  = NaN;
        if ~isnan(B_map(ig))
            if j>=1 && j<=size(Broadcast_GIVD,2)
                truth_val = Broadcast_GIVD(B_map(ig), j);
            end
        end

        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);
        cos_lat1 = cos(lat1*DEG2RAD); % 预先计算，不影响数值

        % 进度提示（保持在外层 epoch 频率，不影响结果）
        if ~isempty(dq) && mod(ie, 500) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('STK线程 %d | epoch %d/%d ', tid, ie, num_epochs));
        end

        % ====================== ST 分支（保持原逻辑） ======================
        lat_cand = pf_lat(time_idx_win);
        lon_cand = pf_lon(time_idx_win);

        % 大圆距离（向量化实现，但与原表达等价）
        dlat = (lat_cand - lat1)*DEG2RAD;
        dlon = (lon_cand - lon1)*DEG2RAD;
        a = sin(0.5*dlat).^2 + cos_lat1 .* cos(lat_cand*DEG2RAD) .* sin(0.5*dlon).^2;
        Dist_all = R_earth .* (2*atan2(sqrt(a), sqrt(1-a)));

        selected = (Dist_all <= R_max & Dist_all > 0);
        idsearch_initial = time_idx_win(selected);
        if numel(idsearch_initial) < N_min
            continue;
        end

        % Step1：ST变差函数排序（保持 sort 升序与原代码完全一致）
        s     = Dist_all(selected);
        tdiff = abs(pf_epoch(idsearch_initial) - est_time);
        vvals = stpacked(ST_ST2model, s, tdiff);
        [sorted_variogram, ord] = sort(vvals, 'ascend'); % 与原一致
        idsearch_initial = idsearch_initial(ord);
        variogram_values = sorted_variogram;

        % 计算极值与阈值（保持原逻辑）
        max_variogram = max(variogram_values);
        min_variogram = min(variogram_values);
        PF_samples    = PF(idsearch_initial, :);
        target_point  = [est_time, ig, lat1, lon1, truth_val]; % 注意：用 j（与原一致）

        if numel(idsearch_initial) >= N_target
            variogram_des = variogram_values(N_target);
        else
            variogram_des = max_variogram;
        end

        % Mode 1: Fixed width（保持原循环顺序与步进）
        R_in_values = linspace(min_variogram-1e-10, variogram_des+1e-10, Rin_numbers+1);
        width_all   = (0.1:0.2:0.9) .* (max_variogram - min_variogram);
        for r = 1:length(R_in_values)-1
            R_in_left = R_in_values(r);
            for Rw = width_all
                R_in  = R_in_left;
                R_out = R_in_left + Rw;
                result_once = STdepriveAndEstimate2(R_min, R_max, N, R_in, R_out, ...
                                PF_samples, target_point, ST_ST2model, wgs84, height, R_earth, 'fixed width');
                if isnan(result_once), continue; end
                save_row = save_row + 1;
                results(save_row, :) = [target_point result_once];
            end
        end

        % Mode 3: Temporal delete（保持原循环顺序与步进）
        R_in_values2 = 0:30:300;
        Rw_values2   = (0.1:0.2:0.9) * 300;
        for r = 1:length(R_in_values2)
            R_in = R_in_values2(r);
            for w = 1:numel(Rw_values2)
                R_out = R_in + Rw_values2(w);
                result_once = STdepriveAndEstimate2(R_min, R_max, N, R_in, R_out, ...
                                PF_samples, target_point, ST_ST2model, wgs84, height, R_earth, 'temporal delete');
                if isnan(result_once), continue; end
                save_row = save_row + 1;
                results(save_row, :) = [target_point result_once];
            end
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
