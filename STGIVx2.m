function [results, CONDs] = STGIVx2(R, N, IPPVD, samp_int, ST_STmodels, type, IGPs, Broadcast_GIVD, Broadcast_IGPs, dq, monothreat)
% STGIVx — Spatio-Temporal estimation on fixed IGP grid points, with Broadcast GIVD as "truth"
%
% 正常：若该历元存在 STmodel2 → 尝试 STkriging；
% 失败或不可用：统一走 UnK_threat 兜底（仅用当前历元）。
%
% outputs（每行17列）:
% [ epoch, Broadcast_GIVD_truth, truth_val2, IGP_lat, IGP_lon, UnK_result(1:6), ...
%   Rfit, RCM, R_T, TCM, sigma_undersampled, sigma_GIVE ]

% ----------- 参数与默认值 -----------
if nargin < 6 || isempty(type),      type = 'GIVD'; end
if nargin < 7 || isempty(IGPs),      error('IGPs (M×2) is required: [lat, lon].'); end
if nargin < 8 || isempty(Broadcast_GIVD) || nargin < 9 || isempty(Broadcast_IGPs)
    error('Broadcast_GIVD (M_b×2880) and Broadcast_IGPs (M_b×2) are required.');
end
if nargin < 10, dq = []; end
if nargin < 11, monothreat = []; end

R_min = R(1); %#ok<NASGU> % 仅保留（未用到时也不报错）
R_max = R(2);
N_min = N(1);
N_target = N(2);

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350000/1000; % km
R_earth = 6378.137;    % km

% === Default UnK 模型（用于 UnK_threat 或 ST 缺省场景的排序距离核）===
defaultUnK.func        = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
defaultUnK.sill        = 1;           % total^2
defaultUnK.model       = 'exponential';
defaultUnK.nugget      = 0.3 * 0.3;   % nom^2
defaultUnK.partialsill = 1 - 0.3 * 0.3;
defaultUnK.range       = 8000;        % km （空间相关距离，用于排序时的尺度）
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
CONDs   = zeros(num_epochs * num_igps, 3);
save_row = 0;

% GIVE 模式分箱
Rfit_edges = 0:10:400;    % km
TCM_edges  = 0:0.025:1;
RCM_edges  = 0:0.025:1;

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

    % 该历元的 STmodel2（若不存在或为空 → 后续将 fallback）
    model_idx = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    has_ST = (~isempty(model_idx) && ~isempty(ST_STmodels(model_idx).model));
    if has_ST
        ST_ST2model = ST_STmodels(model_idx).model;
    else
        ST_ST2model = []; % 显式为空用于判断
    end

    % 时间窗（ST 用过去 300s；UnK_fallback 仅用当前历元）
    time_range   = est_time - 300 : est_time;
    time_idx_win = find(ismember(IPPVD(:,1), time_range));  % 窗口索引
    time_idx_exact = find(IPPVD(:,1) == est_time);          % 当前历元索引

    % 日内 30s 索引 j（列1=00:00:00）
    j  = epoch_to_30s_index(est_time);
    jj = epoch_to_30s_index(est_time) - 1; % 你原代码的第二个参考列
    if jj < 1, jj = 1; end

    for ig = 1:num_igps
        % Broadcast GIVD 真值
        truth_val  = NaN;
        truth_val2 = NaN;
        if ~isnan(j)
            m = B_map(ig);
            if ~isnan(m)
                truth_val  = Broadcast_GIVD(m, j);
                truth_val2 = Broadcast_GIVD(m, jj);
            end
        end

        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);

        % 进度提示
        if ~isempty(dq) && mod(ig, 50) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('线程 %d | epoch %d/%d | IGP %d/%d', tid, ie, num_epochs, ig, num_igps));
        end

        % —— 先尝试 ST（若可用且窗口非空）——
        pack = [];   % 结构体：{UnK_result, COND, Rfitt, RCM1, R_T, TCM, sigma_undersampled, sigma_GIVE}
        use_ST = has_ST && ~isempty(time_idx_win);

        if use_ST
            % 候选 IPP 的经纬度
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

            if numel(idsearch_initial) >= N_min
                % Step1：按 ST 变差函数排序
                s     = Dist_all(selected);
                tdiff = abs(IPPVD(idsearch_initial,1) - est_time);
                vvals = ST_ST2model.func(s, tdiff);
                [~, ord] = sort(vvals, 'ascend');
                idsearch_initial = idsearch_initial(ord);

                % 截到 N_target（或用全体）
                if numel(idsearch_initial) > N_target
                    use_idx = idsearch_initial(1:N_target);
                else
                    use_idx = idsearch_initial;
                end
                PF_samples = IPPVD(use_idx, :);
            else
                PF_samples = [];
            end

            if ~isempty(PF_samples) && size(PF_samples,1) >= N_min
                % IGP—样本时空量
                dlat = deg2rad(PF_samples(:,3) - lat1);
                dlon = deg2rad(PF_samples(:,4) - lon1);
                a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(PF_samples(:,3))).*sin(dlon/2).^2;
                c = 2 * atan2(sqrt(a), sqrt(1-a));
                Dist        = R_earth * c;
                t_remaining = est_time - PF_samples(:,1);

                % Rfitt 规则
                num_cand_used = size(PF_samples,1);
                if num_cand_used >= N_target
                    Rfitt = max(max(Dist), R_min);
                else
                    Rfitt = R_max;
                end
                R_T = max(t_remaining);

                % ENU
                East2  = zeros(size(PF_samples,1),1);
                North2 = zeros(size(PF_samples,1),1);
                for k = 1:size(PF_samples,1)
                    [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                                               lat1, lon1, height, wgs84, 'degrees');
                    East2(k)  = E;  North2(k) = Nn;
                end

                % 距离/时间差矩阵
                num_samples = size(PF_samples,1);
                dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
                [Mgrid, Ngrid] = meshgrid(1:num_samples);
                lat_m = PF_samples(Mgrid(:),3); lon_m = PF_samples(Mgrid(:),4);
                lat_n = PF_samples(Ngrid(:),3); lon_n = PF_samples(Ngrid(:),4);
                dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n, R_earth), num_samples, num_samples);
                t_m  = PF_samples(Mgrid(:),1); t_n  = PF_samples(Ngrid(:),1);
                t_nn = reshape(abs(t_m - t_n), num_samples, num_samples);

                % RCM/TCM
                lati_mean   = mean(PF_samples(:,3));
                long_mean   = mean(PF_samples(:,4));
                Rcentroid1  = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
                RCM1        = Rcentroid1 / max(Rfitt, eps);
                if R_T == 0, TCM = 0; else, TCM = mean(t_remaining) / R_T; end

                % STkriging（truth_val 仅用于诊断指标）
                data_values = PF_samples(:,5);
                [X0, Y0] = meshgrid(0);
                [UnK_result_ST, COND_ST] = STkriging(ST_ST2model, East2, North2, t_remaining, data_values, ...
                                                     X0, Y0, 0, dist_nn, dist_n0, t_nn, truth_val);

                if ~isempty(UnK_result_ST)
                    % GIVE（3D）
                    sigma_undersampled = NaN; sigma_GIVE = NaN;
                    if strcmpi(type, 'GIVE') && ~isempty(monothreat)
                        Rfit_bin = discretize(Rfitt, Rfit_edges);
                        RCM_bin  = discretize(RCM1,  RCM_edges);
                        TCM_bin  = discretize(TCM,   TCM_edges);
                        if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                            sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                            sigma_GIVE = sqrt(sigma_undersampled.^2 + UnK_result_ST(6).^2);
                        end
                    end

                    % 打包（ST 成功）
                    pack = struct('UnK_result',UnK_result_ST,'COND',COND_ST, ...
                                  'Rfitt',Rfitt,'RCM1',RCM1,'R_T',R_T,'TCM',TCM, ...
                                  'sigma_undersampled',sigma_undersampled,'sigma_GIVE',sigma_GIVE);
                end
            end
        end

        % —— 若 ST 不可用或失败：统一走 UnK fallback ——
        if isempty(pack)
            pack = do_UnK_fallback(est_time, lat1, lon1, IPPVD, time_idx_exact, ...
                                   R_min, R_max, N_min, N_target, ...
                                   wgs84, height, R_earth, ...
                                   defaultUnK, defaultfit, ...
                                   monothreat, type, Rfit_edges, RCM_edges, TCM_edges);
            if isempty(pack), continue; end  % fallback 也失败就跳过
        end

        % —— 写结果（统一一处）——
        save_row = save_row + 1;
        results(save_row, :) = [ ...
            est_time, truth_val, truth_val2, lat1, lon1, ...
            pack.UnK_result(:).', ...
            pack.Rfitt, pack.RCM1, pack.R_T, pack.TCM, ...
            pack.sigma_undersampled, pack.sigma_GIVE ...
        ];
        CONDs(save_row, :) = pack.COND(:).';
    end
end

% 裁剪
results = results(1:save_row, :);
CONDs   = CONDs(1:save_row, :);

% ====================== 子函数：UnK 兜底 ======================
function pack = do_UnK_fallback(est_time, lat1, lon1, IPPVD, time_idx_exact, ...
                                R_min, R_max, N_min, N_target, ...
                                wgs84, height, R_earth, ...
                                defaultUnK, defaultfit, ...
                                monothreat, type, Rfit_edges, RCM_edges, TCM_edges)
    pack = [];
    PF_now = IPPVD(time_idx_exact, :);
    if isempty(PF_now), return; end

    % IGP—所有 IPP 的距离
    lat2 = PF_now(:,3); lon2 = PF_now(:,4);
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    Dist_full = R_earth * c;

    % 排序并去掉可能重合点
    [Dist_sorted, idx_sort] = sort(Dist_full, 'ascend'); %#ok<ASGLU>
    if ~isempty(Dist_sorted), Dist_sorted(1) = []; end

    % 规则确定 Rfitt
    num_IPPs_le_Rmin = sum(Dist_sorted <= R_min);
    if num_IPPs_le_Rmin >= N_target
        Rfitt = R_min;
    elseif numel(Dist_sorted) >= N_target && Dist_sorted(N_target) <= R_max
        Rfitt = Dist_sorted(N_target);
    else
        Rfitt = R_max;
        if sum(Dist_full <= Rfitt) - any(Dist_full==0) < N_min
            return; % Rmax 内样本不足
        end
    end

    % 选取样本
    sel        = (Dist_full <= Rfitt) & (Dist_full > 0);
    PF_samples = PF_now(sel, :);
    if size(PF_samples,1) < N_min, return; end

    % ENU
    East2  = zeros(size(PF_samples,1),1);
    North2 = zeros(size(PF_samples,1),1);
    for k = 1:size(PF_samples,1)
        [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                                   lat1, lon1, height, wgs84, 'degrees');
        East2(k)  = E;  North2(k) = Nn;
    end

    % 距离矩阵（无时间维）
    num_samples = size(PF_samples,1);
    dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
    [Mgrid, Ngrid] = meshgrid(1:num_samples);
    lat_m = PF_samples(Mgrid(:),3); lon_m = PF_samples(Mgrid(:),4);
    lat_n = PF_samples(Ngrid(:),3); lon_n = PF_samples(Ngrid(:),4);
    dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n, R_earth), num_samples, num_samples);

    % RCM/时间指标
    lati_mean   = mean(PF_samples(:,3));
    long_mean   = mean(PF_samples(:,4));
    Rcentroid1  = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
    RCM1        = Rcentroid1 / max(Rfitt, eps);
    R_T = 0; TCM = 0;

    % UnK 估计
    data_values = PF_samples(:,5);
    [X0, Y0]    = meshgrid(0);
    UnK_result  = Unk_threat(defaultUnK, East2, North2, data_values, X0, Y0, defaultfit, dist_nn, dist_n0);

    % GIVE/GIVD（兼容 2D/3D monothreat）
    sigma_undersampled = NaN; sigma_GIVE = NaN;
    if ~isempty(monothreat)
        if ndims(monothreat) == 2
            Rfit_bin = discretize(Rfitt, Rfit_edges);
            RCM_bin  = discretize(RCM1,  RCM_edges);
            if ~isnan(Rfit_bin) && ~isnan(RCM_bin)
                sigma_undersampled = monothreat(RCM_bin, Rfit_bin);
            end
        else
            Rfit_bin = discretize(Rfitt, Rfit_edges);
            RCM_bin  = discretize(RCM1,  RCM_edges);
            TCM_bin  = 1; % 假设第1层为 TCM=0
            if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
            end
        end
        if isfinite(sigma_undersampled)
            sigma_GIVE = sqrt(sigma_undersampled.^2 + UnK_result(6).^2);
        end
    end

    pack = struct('UnK_result',UnK_result,'COND',[NaN,NaN,NaN], ...
                  'Rfitt',Rfitt,'RCM1',RCM1,'R_T',R_T,'TCM',TCM, ...
                  'sigma_undersampled',sigma_undersampled,'sigma_GIVE',sigma_GIVE);
end

end

% === 辅助函数：列1=00:00:00 → j = floor(sod/30)+1 ===
function j = epoch_to_30s_index(epoch_seconds)
    if ~isfinite(epoch_seconds), j = NaN; return; end
    sod = mod(epoch_seconds, 86400);
    if sod < 0, sod = sod + 86400; end
    j = floor(sod / 30) + 1; % 0→1, 30→2, ..., 86370→2880
end
