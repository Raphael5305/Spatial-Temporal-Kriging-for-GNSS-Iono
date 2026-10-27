function [results, CONDs] = UnKGIVx(R, N, IPPVD, samp_int, type, IGPs, Broadcast_GIVD, Broadcast_IGPs, dq, monothreat)
% UnKGIVx — Spatial (UnK) estimation on fixed IGP grid points, with Broadcast GIVD as "truth"
%
% 仅使用空间 UnK（无 STkriging），每个历元独立，用 WAAS 四步规则确定 Rfitt，
% 在 Rfitt 圆内取样本做 Unk_threat 估计，并输出 GIVE/GIVD 相关指标。
%
% outputs（每行17列）:
% [ epoch, Broadcast_GIVD_truth, Broadcast_GIVD_truth_prev, IGP_lat, IGP_lon, ...
%   UnK_result(1:6), Rfit, RCM, R_T(=0), TCM(=0), sigma_undersampled, sigma_GIVE ]

% ----------- 参数与默认值 -----------
if nargin < 5 || isempty(type),      type = 'GIVD'; end
if nargin < 6 || isempty(IGPs),      error('IGPs (M×2) is required: [lat, lon].'); end
if nargin < 7 || isempty(Broadcast_GIVD) || nargin < 8 || isempty(Broadcast_IGPs)
    error('Broadcast_GIVD (M_b×2880) and Broadcast_IGPs (M_b×2) are required.');
end
if nargin < 9,  dq = []; end
if nargin < 10, monothreat = []; end

R_min = R(1);
R_max = R(2);
N_min = N(1);
N_target = N(2);

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350000/1000; % km
R_earth = 6378.137;    % km

% === Default UnK 模型 ===
defaultUnK.func        = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
defaultUnK.sill        = 1;           % total^2
defaultUnK.model       = 'exponential';
defaultUnK.nugget      = 0.3 * 0.3;   % nom^2
defaultUnK.partialsill = 1 - 0.3 * 0.3;
defaultUnK.range       = 8000;        % km
defaultfit             = [defaultUnK.range, defaultUnK.partialsill, defaultUnK.nugget];

% --- 检查 Broadcast 输入 ---
if size(Broadcast_IGPs,2) ~= 2
    error('Broadcast_IGPs must be M_b×2 (lat, lon).');
end
if size(Broadcast_GIVD,2) ~= 2880
    error('Broadcast_GIVD must have 2880 columns (one day at 30s cadence).');
end
if size(Broadcast_GIVD,1) ~= size(Broadcast_IGPs,1)
    error('Rows of Broadcast_GIVD must match rows of Broadcast_IGPs.');
end

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

% 预分配
results = zeros(num_epochs * num_igps, 17);
CONDs   = zeros(num_epochs * num_igps, 3);  % UnK 无条件数输出，置 NaN
save_row = 0;

% GIVE 分箱
Rfit_edges = 0:10:400;    % km
RCM_edges  = 0:0.025:1;

% --- IGP → Broadcast_IGPs 映射（精确匹配） ---
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
for ie = 1:num_epochs
    est_time = est_epochs(ie);

    % 当前历元索引
    time_idx_exact = find(IPPVD(:,1) == est_time);
    if isempty(time_idx_exact)
        continue;
    end

    % 日内 30s 索引 j（列1=00:00:00）
    j  = epoch_to_30s_index(est_time);
    jj = j - 1; if jj < 1, jj = NaN; end  % 上一个 30s 桢的索引（若不足一列则 NaN）

    for ig = 1:num_igps
        % Broadcast GIVD 真值
        truth_val = NaN; truth_val2 = NaN;
        m = B_map(ig);
        if ~isnan(m) && ~isnan(j)
            truth_val  = Broadcast_GIVD(m, j);
            if ~isnan(jj), truth_val2 = Broadcast_GIVD(m, jj); end
        end

        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);

        % 进度提示
        if ~isempty(dq) && mod(ig, 50) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('线程 %d | epoch %d/%d | IGP %d/%d', tid, ie, num_epochs, ig, num_igps));
        end

        % =============== UnK（当前历元） ===============
        PF_now = IPPVD(time_idx_exact, :);
        if isempty(PF_now), continue; end

        % IGP-IPP 大圆距离
        lat2 = PF_now(:,3);
        lon2 = PF_now(:,4);
        dlat = deg2rad(lat2 - lat1);
        dlon = deg2rad(lon2 - lon1);
        a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        Dist_full = R_earth * c;

        % 排序并跳过自身（若存在）
        [Dist_sorted, ord] = sort(Dist_full, 'ascend');
        if ~isempty(Dist_sorted) && Dist_sorted(1) == 0
            Dist_sorted(1) = [];
            ord(1) = [];
        end

        % WAAS 四步规则确定 Rfitt
        nRmin = sum(Dist_sorted <= R_min);
        if nRmin >= N_target
            Rfitt = R_min;
        elseif numel(Dist_sorted) >= N_target && Dist_sorted(N_target) <= R_max
            Rfitt = Dist_sorted(N_target);
        else
            Rfitt = R_max;
            if sum(Dist_sorted <= Rfitt) < N_min
                continue; % Rmax 内样本不足
            end
        end

        % 取 Rfitt 内的样本
        sel = (Dist_full <= Rfitt) & (Dist_full > 0);
        PF_samples = PF_now(sel, :);
        Dist       = Dist_full(sel);
        if size(PF_samples,1) < N_min, continue; end

        % ENU
        East2  = zeros(size(PF_samples,1),1);
        North2 = zeros(size(PF_samples,1),1);
        for k = 1:size(PF_samples,1)
            [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                                       lat1, lon1, height, wgs84, 'degrees');
            East2(k)  = E;
            North2(k) = Nn;
        end

        % 距离矩阵（无时间维）
        num_samples = size(PF_samples,1);
        dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
        [Mgrid, Ngrid] = meshgrid(1:num_samples);
        lat_m = PF_samples(Mgrid(:),3);
        lon_m = PF_samples(Mgrid(:),4);
        lat_n = PF_samples(Ngrid(:),3);
        lon_n = PF_samples(Ngrid(:),4);
        dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n, R_earth), num_samples, num_samples);

        % RCM（无时间维）
        lati_mean = mean(PF_samples(:,3));
        long_mean = mean(PF_samples(:,4));
        Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
        RCM1 = Rcentroid1 / max(Rfitt, eps);

        % UnK 估计
        data_values = PF_samples(:,5);
        [X0, Y0] = meshgrid(0);
        UnK_result = Unk_threat(defaultUnK, East2, North2, data_values, X0, Y0, defaultfit, dist_nn, dist_n0);
        % UnK_result: [zi, chi, chi_irreg, R_irreg, sigma_IPP, sigma_inflated]

        % GIVE/GIVD
        sigma_undersampled = NaN; sigma_GIVE = NaN;
        if ~isempty(monothreat)
            if ndims(monothreat) == 2
                % monothreat(RCM_bin, Rfit_bin)
                Rfit_bin = discretize(Rfitt, Rfit_edges);
                RCM_bin  = discretize(RCM1, RCM_edges);
                if ~isnan(Rfit_bin) && ~isnan(RCM_bin)
                    sigma_undersampled = monothreat(RCM_bin, Rfit_bin);
                end
            else
                % 若提供三维表，当前无时间维，保守取 TCM=0 的切片
                Rfit_bin = discretize(Rfitt, Rfit_edges);
                RCM_bin  = discretize(RCM1, RCM_edges);
                TCM_bin  = 1;
                if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                    sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                end
            end
            if isfinite(sigma_undersampled)
                sigma_GIVE = sqrt(sigma_undersampled.^2 + UnK_result(6).^2);
            end
        end

        % 本方法无时间维：R_T = 0, TCM = 0；COND 置 NaN
        R_T = 0; TCM = 0; COND = [NaN, NaN, NaN];

        % 写结果
        save_row = save_row + 1;
        results(save_row, :) = [ ...
            est_time, truth_val, truth_val2, lat1, lon1, ...
            UnK_result(:).', ...
            Rfitt, RCM1, R_T, TCM, ...
            sigma_undersampled, sigma_GIVE ...
        ];
        CONDs(save_row, :) = COND(:).';
    end
end

% 裁剪
results = results(1:save_row, :);
CONDs   = CONDs(1:save_row, :);
end

% === 辅助函数：列1=00:00:00 → j = floor(sod/30)+
function j = epoch_to_30s_index(epoch_seconds)
    if ~isfinite(epoch_seconds), j = NaN; return; end
    sod = mod(epoch_seconds, 86400);
    if sod < 0, sod = sod + 86400; end
    j = floor(sod / 30) + 1; % 0→1, 30→2, ..., 86370→2880
end
