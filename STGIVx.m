function [results, CONDs] = STGIVx(R, N, IPPVD, samp_int, ST_STmodels, type, IGPs, Broadcast_GIVD, Broadcast_IGPs, dq, monothreat)
% STGIVx — Spatio-Temporal estimation on fixed IGP grid points, with Broadcast GIVD as "truth"
%
% 正常：若该历元存在 STmodel2 → 用 STkriging；
% 兜底：若该历元 STmodel2 为空 → 用 UnK_threat（与 UnKGIVE 思路一致），极少数情况
%
% outputs（每行16列）:
% [ epoch, Broadcast_GIVD_truth, IGP_lat, IGP_lon, UnK_result(1:6), ...
%   Rfit, RCM, R_T, TCM, sigma_undersampled, sigma_GIVE ]

% ----------- 参数与默认值 -----------
if nargin < 6 || isempty(type),      type = 'GIVD'; end
if nargin < 7 || isempty(IGPs),      error('IGPs (M×2) is required: [lat, lon].'); end
if nargin < 8 || isempty(Broadcast_GIVD) || nargin < 9 || isempty(Broadcast_IGPs)
    error('Broadcast_GIVD (M_b×2880) and Broadcast_IGPs (M_b×2) are required.');
end
if nargin < 10, dq = []; end
if nargin < 11, monothreat = []; end

R_min = R(1); %#ok<NASGU> % 仅保留
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
for ie = 1:num_epochs
    est_time = est_epochs(ie);

    % 该历元的 STmodel2（若不存在或为空 → 走 UnK_threat 分支）
    model_idx = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    model_idxave=ceil(model_idx/10)*10;
    model_idxave=min(2870,model_idxave);
    has_ST = (~isempty(model_idx) && ~isempty(ST_STmodels(model_idx).packed));
    has_average = (~isempty(model_idx) && ~isempty(ST_STmodels(model_idxave).packavegrage));
    ST_model     = [];
    ST_model_avg = [];
    if has_ST,      ST_model     = ST_STmodels(model_idx).packed;             end
    if has_average, ST_model_avg = ST_STmodels(model_idxave).packavegrage;    end

    % 时间窗候选（ST 分支用过去300s；UnK_threat 将改为仅用当前历元）
    time_range = est_time - 300 : est_time;
    time_idx_win   = find(ismember(IPPVD(:,1), time_range));  % 窗口索引
    time_idx_exact = find(IPPVD(:,1) == est_time);            % 当前历元索引

    if isempty(time_idx_win) 
        % ST需要窗口；若窗口空，直接跳过这个历元
        continue;
    end

    if isempty(time_idx_exact) && ~(has_ST || has_average), continue; end %UNKm模式但无当前历元


    % 日内 30s 索引 j（列1=00:00:00）
    j = epoch_to_30s_index(est_time);
    jj=epoch_to_30s_index(est_time)-1;%判断哪个时间是对的

    for ig = 1:num_igps
        % Broadcast GIVD 真值（第2列 & 传给 STkriging 用于诊断）
        truth_val = NaN;
        truth_val2=NaN;
        if ~isnan(j)
            m = B_map(ig);
            if ~isnan(m), truth_val = Broadcast_GIVD(m, j); truth_val2 = Broadcast_GIVD(m, jj);end
        end

        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);

        % 进度提示
        if ~isempty(dq) && mod(ig, 50) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('线程 %d | epoch %d/%d | IGP %d/%d', tid, ie, num_epochs, ig, num_igps));
        end

        if has_ST|| has_average
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
            id_all = time_idx_win(selected);
            if numel(id_all) < N_min %样本不足直接跳过
                continue;
            end

            % Step1：以 ST 变差函数排序
            s = Dist_all(selected);
            tdiff = abs(IPPVD(id_all,1) - est_time);

            % 最佳候选占位
            got_best = false;
            best_err = Inf;          % 只用 |zi - truth_val|
            best_vec = [];           % UnK_result(1:6)
            best_COND = [NaN,NaN,NaN];
            best_Rfitt = NaN; best_R_T = NaN; best_TCM = NaN; best_RCM1 = NaN;

         % ===== 策略A：用 ST 模型排序选样（若 has_ST）=====
         if has_ST
             v_st = stpacked(ST_model, s, tdiff);
                [~, ord_st] = sort(v_st, 'ascend');
            id_st = id_all(ord_st);

            % 截取前 N_target
            if numel(id_st) > N_target
                use_idx = id_st(1:N_target);
            else%使用全部IPPs：Nmin-Ntarget
                use_idx = id_st;
            end
            PF_samples = IPPVD(use_idx, :);
            if size(PF_samples,1) < N_min, continue; end


            % IGP 与样本的距离/时间差
            dlat = deg2rad(PF_samples(:,3) - lat1);
            dlon = deg2rad(PF_samples(:,4) - lon1);
            a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(PF_samples(:,3))).*sin(dlon/2).^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            Dist = R_earth * c; % km
            t_remaining = est_time - PF_samples(:,1);

            % ====== 新：Rfitt 由“已选样本”的几何外延决定（并不低于 R_min；欠采样则置 R_max）======
            num_cand_used = numel(id_st);  % ≤ 所有候选(≤R_max)；>=N_min已保证
            if num_cand_used >= N_target
                % 正常：选满 N_target，则 Rfitt = max(选中样本距离, R_min)
                Rfitt_A = max(max(Dist), R_min);
            else
                % 欠采样：R_max 内总样本数在 [N_min, N_target) ⇒ Rfitt = R_max（WAAS 规则）
                Rfitt_A = R_max;
            end
            R_T_A   = max(t_remaining);

            % East/North
            East2  = zeros(size(PF_samples,1),1);
            North2 = zeros(size(PF_samples,1),1);
            for k = 1:size(PF_samples,1)
                [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                                           lat1, lon1, height, wgs84, 'degrees');
                East2(k)  = E;
                North2(k) = Nn;
            end

            % 距离矩阵与时间差矩阵
            num_samples = size(PF_samples,1);
            dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
            [Mgrid, Ngrid] = meshgrid(1:num_samples);
            lat_m = PF_samples(Mgrid(:),3);
            lon_m = PF_samples(Mgrid(:),4);
            lat_n = PF_samples(Ngrid(:),3);
            lon_n = PF_samples(Ngrid(:),4);
            dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n, R_earth), num_samples, num_samples);
            t_m  = PF_samples(Mgrid(:),1);
            t_n  = PF_samples(Ngrid(:),1);
            t_nn = reshape(abs(t_m - t_n), num_samples, num_samples);

            % RCM / TCM
            lati_mean = mean(PF_samples(:,3));
            long_mean = mean(PF_samples(:,4));
            Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
            RCM1_A = Rcentroid1 / max(Rfitt_A, eps);
            if R_T_A == 0, TCM_A = 0; else, TCM_A = mean(t_remaining) / R_T_A; end

            % STkriging（用 truth_val 仅作诊断，不参与权重）
            data_values = PF_samples(:,5);
            [X0, Y0] = meshgrid(0);
            target_truth = truth_val;
            target_truth2 = truth_val2;
            [uk_A1, COND_A1] = STkriging(ST_model, East2, North2, t_remaining, data_values, ...
                                           X0, Y0, 0, dist_nn, dist_n0, t_nn, target_truth);
            if ~isempty(uk_A1)
                cur_err = abs(uk_A1(1) - truth_val);
                if isfinite(cur_err)
                    if cur_err < best_err
                        got_best = true; best_err = cur_err;
                        best_vec = uk_A1; best_COND = COND_A1;
                        best_Rfitt = Rfitt_A; best_R_T = R_T_A; best_TCM = TCM_A; best_RCM1 = RCM1_A;
                    end
                elseif ~got_best
                    got_best = true; best_err = cur_err;
                    best_vec = uk_A1; best_COND = COND_A1;
                    best_Rfitt = Rfitt_A; best_R_T = R_T_A; best_TCM = TCM_A; best_RCM1 = RCM1_A;
                end
            end
            % A2) ST选样 + average模型（若 has_average）
            if has_average
                [uk_A2, COND_A2] = STkriging(ST_model_avg, East2, North2, t_remaining, data_values, ...
                    X0, Y0, 0, dist_nn, dist_n0, t_nn, truth_val);
                if ~isempty(uk_A2)
                    cur_err = abs(uk_A2(1) - truth_val);
                    if isfinite(cur_err)
                        if cur_err < best_err
                            got_best = true; best_err = cur_err;
                            best_vec = uk_A2; best_COND = COND_A2;
                            best_Rfitt = Rfitt_A; best_R_T = R_T_A; best_TCM = TCM_A; best_RCM1 = RCM1_A;
                        end
                    elseif ~got_best
                        got_best = true; best_err = cur_err;
                        best_vec = uk_A2; best_COND = COND_A2;
                        best_Rfitt = Rfitt_A; best_R_T = R_T_A; best_TCM = TCM_A; best_RCM1 = RCM1_A;
                    end
                end
            end
         end


% ===== 策略B：用 average 模型排序选样 + average 模型（若 has_average）=====
            if has_average
                v_av = stpacked(ST_model_avg, s, tdiff);
                [~, ord_av] = sort(v_av, 'ascend');
                id_av_sel = id_all(ord_av);
                if numel(id_av_sel) > N_target, id_av_sel = id_av_sel(1:N_target); end
                PF_av = IPPVD(id_av_sel, :);

                % 选样后不足：直接跳过
                if size(PF_av,1) < N_min,continue;end
                    dlat = deg2rad(PF_av(:,3) - lat1);
                    dlon = deg2rad(PF_av(:,4) - lon1);
                    a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(PF_av(:,3))).*sin(dlon/2).^2;
                    c = 2 * atan2(sqrt(a), sqrt(1-a));
                    Dist = R_earth * c;
                    t_remaining = est_time - PF_av(:,1);

                    if numel(id_av_sel) >= N_target
                        Rfitt_B = max(max(Dist), R_min);
                    else
                        Rfitt_B = R_max;
                    end
                    R_T_B = max(t_remaining);

                    East2  = zeros(size(PF_av,1),1);
                    North2 = zeros(size(PF_av,1),1);
                    for kk = 1:size(PF_av,1)
                        [E, Nn, ~] = geodetic2enu(PF_av(kk,3), PF_av(kk,4), height, ...
                                                   lat1, lon1, height, wgs84, 'degrees');
                        East2(kk)  = E;  North2(kk) = Nn;
                    end

                    ns = size(PF_av,1);
                    dist_n0 = distance('gc', PF_av(:,3), PF_av(:,4), lat1, lon1, R_earth);
                    [Mg, Ng] = meshgrid(1:ns);
                    dist_nn = reshape(distance('gc', PF_av(Mg(:),3), PF_av(Mg(:),4), ...
                                                    PF_av(Ng(:),3), PF_av(Ng(:),4), R_earth), ns, ns);
                    t_nn = reshape(abs(PF_av(Mg(:),1) - PF_av(Ng(:),1)), ns, ns);

                    lati_mean = mean(PF_av(:,3)); long_mean = mean(PF_av(:,4));
                    Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
                    RCM1_B = Rcentroid1 / max(Rfitt_B, eps);
                    if R_T_B == 0, TCM_B = 0; else, TCM_B = mean(t_remaining) / R_T_B; end

                    data_values = PF_av(:,5);
                    [X0, Y0] = meshgrid(0);

                    [uk_B, COND_B] = STkriging(ST_model_avg, East2, North2, t_remaining, data_values, ...
                                               X0, Y0, 0, dist_nn, dist_n0, t_nn, truth_val);
                    if ~isempty(uk_B)
                        cur_err = abs(uk_B(1) - truth_val);
                        if isfinite(cur_err)
                            if cur_err < best_err
                                got_best = true; best_err = cur_err;
                                best_vec = uk_B; best_COND = COND_B;
                                best_Rfitt = Rfitt_B; best_R_T = R_T_B; best_TCM = TCM_B; best_RCM1 = RCM1_B;
                            end
                        elseif ~got_best
                            got_best = true; best_err = cur_err;
                            best_vec = uk_B; best_COND = COND_B;
                            best_Rfitt = Rfitt_B; best_R_T = R_T_B; best_TCM = TCM_B; best_RCM1 = RCM1_B;
                        end
                    end

            end
 
            % 只在选出最佳候选后，计算一次 sigma_undersampled / sigma_GIVE 并写结果
            if got_best
                sigma_undersampled = NaN; sigma_GIVE = NaN;
                if strcmpi(type,'GIVE') && ~isempty(monothreat)
                    Rfit_bin = discretize(best_Rfitt, Rfit_edges);
                    RCM_bin  = discretize(best_RCM1, RCM_edges);
                    TCM_bin  = discretize(best_TCM,  TCM_edges);
                    if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                        sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                        if numel(best_vec) >= 6 && isfinite(best_vec(6))
                            % 按你原 ST 分支写法（不对第6列再平方）
                            sigma_GIVE = sqrt(sigma_undersampled.^2 + best_vec(6));
                        end
                    end
                end

                save_row = save_row + 1;
                results(save_row, :) = [ ...
                    est_time, truth_val, truth_val2, lat1, lon1, ...
                    best_vec(:).', ...
                    best_Rfitt, best_RCM1, best_R_T, best_TCM, ...
                    sigma_undersampled, sigma_GIVE ...
                ];
                CONDs(save_row, :) = best_COND(:).';
                continue; % 本 IGP 已完成，跳过 UnK
            else
                % 没有任何有效候选：直接跳过（不转入 UnK）
                continue;
            end
        else
            % ====================== UnK_threat 分支（STmodel2 为空时） ======================
            % 仅用“当前历元”的 IPP（与 UnKGIVE 一致）
            PF_now = IPPVD(time_idx_exact, :);
            if isempty(PF_now), continue; end

            % 计算所有 IPP 与 IGP 的距离
            lat2 = PF_now(:,3);
            lon2 = PF_now(:,4);
            dlat = deg2rad(lat2 - lat1);
            dlon = deg2rad(lon2 - lon1);
            a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            Dist_full = R_earth * c;

            % 排序并跳过自身（距离最小那个可能是 IGP==IPP 的点）
            [Dist_sorted, idx_sort] = sort(Dist_full, 'ascend');
            if ~isempty(Dist_sorted), Dist_sorted(1) = []; idx_sort(1) = []; end

            % 确定 R_fit，参考"SBAS performance improvement with a new undersampled ionosphere threat model based on relative coverage metric"
            num_IPPs_le_Rmin = sum(Dist_sorted <= R_min);
            if num_IPPs_le_Rmin >= N_target
                 % ① Rmin 内样本数已达 N_target
                Rfitt = R_min;
            elseif numel(Dist_sorted) >= N_target && Dist_sorted(N_target) <= R_max
                % ② 在不超过 Rmax 的范围内，可以按距离依次加入直到凑满 N_target
                %    R_fit 取“第 N_target 个被纳入的点”的距离（<= R_max）
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

            % East/North
            East2  = zeros(size(PF_samples,1),1);
            North2 = zeros(size(PF_samples,1),1);
            for k = 1:size(PF_samples,1)
                [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                                           lat1, lon1, height, wgs84, 'degrees');
                East2(k)  = E;
                North2(k) = Nn;
            end

            % 距离矩阵（当前历元 → 没有时间维度）
            num_samples = size(PF_samples,1);
            dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
            [Mgrid, Ngrid] = meshgrid(1:num_samples);
            lat_m = PF_samples(Mgrid(:),3);
            lon_m = PF_samples(Mgrid(:),4);
            lat_n = PF_samples(Ngrid(:),3);
            lon_n = PF_samples(Ngrid(:),4);
            dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n, R_earth), num_samples, num_samples);

            % RCM（无 T 维）
            lati_mean = mean(PF_samples(:,3));
            long_mean = mean(PF_samples(:,4));
            Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
            RCM1 = Rcentroid1 / max(Rfitt, eps);

            % UnK_threat 估计（需你工程中存在 Unk_threat.m）
            data_values = PF_samples(:,5);
            [X0, Y0] = meshgrid(0);
            % 这里调用与你 UnKGIVE 一致的接口：
            % Unk_threat(vstruct_UnK, East, North, data_values, X, Y, UnK_fit, dist_nn, dist_n0)
            UnK_result = Unk_threat(defaultUnK, East2, North2, data_values, X0, Y0, defaultfit, dist_nn, dist_n0);
            % UnK_result: [zi, chi, chi_irreg, R_irreg, sigma_IPP, sigma_inflated]

            % 兜底分支没有时间维度 → monothreat 可能是二维 (RCM, Rfit) 或三维
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
                    % 三维表：没有时间维，保守取 TCM=0 对应的切片
                    Rfit_bin = discretize(Rfitt, Rfit_edges);
                    RCM_bin  = discretize(RCM1, RCM_edges);
                    TCM      = 0; %#ok<NASGU> % 仅为一致性
                    TCM_bin  = 1; % 假定 monothreat 第一层为 TCM=0
                    if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                        sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                    end
                end
                if isfinite(sigma_undersampled)
                    sigma_GIVE = sqrt(sigma_undersampled.^2 + UnK_result(6).^2);
                end
            end

            % 兜底分支中 R_T/TCM 没有意义（仅该历元），统一置 0
            R_T = 0; TCM = 0;
            COND = [NaN, NaN, NaN]; % 无 STkriging 条件量
        end

        % 写结果
        save_row = save_row + 1;
        results(save_row, :) = [ ...
            est_time, truth_val, truth_val2,lat1, lon1, ...
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

% === 辅助函数：列1=00:00:00 → j = floor(sod/30)+1 ===
function j = epoch_to_30s_index(epoch_seconds)
    if ~isfinite(epoch_seconds), j = NaN; return; end
    sod = mod(epoch_seconds, 86400);
    if sod < 0, sod = sod + 86400; end
    j = floor(sod / 30) + 1; % 0→1, 30→2, ..., 86370→2880
end
