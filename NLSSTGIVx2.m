function [results, CONDs] = NLSSTGIVx2( ...
    R, N, IPPVD, samp_int, ST_STmodels, ST_STmodels2, ST_STmodels3, ...
    type, IGPs, Broadcast_GIVD, Broadcast_IGPs, dq, monothreat)
% NLSSTGIVx — Spatio-Temporal estimation on fixed IGP grid points,
% with Broadcast GIVD as "truth"
%
% 六模型独立策略：
%   - 对 ST1/avg1/ST2/avg2/ST3/avg3 六个模型，各自用 stpacked(s,tdiff) 排序→选 N_target 样本→STkriging；
%   - 六个结果与真值 |zi - truth| 比较，取误差最小者为 best；
%   - 仅当六个模型在该历元均不存在时，才启用 UnK 兜底模式；
%   - 其它逻辑（分箱、sigma_undersampled/GIVE 计算、结果列位次）保持一致。
%
% 输出（每行 17 列）:
% [ epoch, Broadcast_GIVD_truth, Broadcast_GIVD_truth_prev, IGP_lat, IGP_lon, uk_best(1:6)或UnK_result(1:6), ...
%   Rfit, RCM, R_T, TCM, sigma_undersampled, sigma_GIVE ]

% ----------- 参数与默认值 -----------
if nargin < 8 || isempty(type),      type = 'GIVD'; end
if nargin < 9 || isempty(IGPs),      error('IGPs (M×2) is required: [lat, lon].'); end
if nargin < 10 || isempty(Broadcast_GIVD) || nargin < 11 || isempty(Broadcast_IGPs)
    error('Broadcast_GIVD (M_b×2880) and Broadcast_IGPs (M_b×2) are required.');
end
if nargin < 12, dq = []; end
if nargin < 13, monothreat = []; end

R_min = R(1);
R_max = R(2);
N_min = N(1);
N_target = N(2);

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350000/1000; % km
R_earth = 6378.137;    % km

% === Default UnK 模型（兜底用）===
defaultUnK.func        = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
defaultUnK.sill        = 1;
defaultUnK.model       = 'exponential';
defaultUnK.nugget      = 0.3 * 0.3;
defaultUnK.partialsill = 1 - 0.3 * 0.3;
defaultUnK.range       = 8000;                 % km
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

% GIVE 分箱
Rfit_edges = 0:10:400;    % km
TCM_edges  = 0:0.025:1;
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

    % === 第一套 ===
    model_idx    = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    model_idxave = []; if ~isempty(model_idx), model_idxave = ceil(model_idx/10)*10; end
    if ~isempty(model_idxave), model_idxave = min(2870,model_idxave); end
    has_ST       = (~isempty(model_idx)    && ~isempty(ST_STmodels(model_idx).packed));
    has_average  = (~isempty(model_idxave) && ~isempty(ST_STmodels(model_idxave).packavegrage));
    ST_model     = []; if has_ST,      ST_model     = ST_STmodels(model_idx).packed;            end
    ST_model_avg = []; if has_average, ST_model_avg = ST_STmodels(model_idxave).packavegrage;   end

    % === 第二套 ===
    model_idx2    = find([ST_STmodels2.epoch_time] == est_time, 1, 'first');
    model_idxave2 = []; if ~isempty(model_idx2), model_idxave2 = ceil(model_idx2/10)*10; end
    if ~isempty(model_idxave2), model_idxave2 = min(2870,model_idxave2); end
    has_ST2       = (~isempty(model_idx2)    && ~isempty(ST_STmodels2(model_idx2).packed));
    has_average2  = (~isempty(model_idxave2) && ~isempty(ST_STmodels2(model_idxave2).packavegrage));
    ST_model2      = []; if has_ST2,      ST_model2      = ST_STmodels2(model_idx2).packed;            end
    ST_model_avg2  = []; if has_average2, ST_model_avg2  = ST_STmodels2(model_idxave2).packavegrage;   end

    % === 第三套（新增）===
    model_idx3    = find([ST_STmodels3.epoch_time] == est_time, 1, 'first');
    model_idxave3 = []; if ~isempty(model_idx3), model_idxave3 = ceil(model_idx3/10)*10; end
    if ~isempty(model_idxave3), model_idxave3 = min(2870,model_idxave3); end
    has_ST3       = (~isempty(model_idx3)    && ~isempty(ST_STmodels3(model_idx3).packed));
    has_average3  = (~isempty(model_idxave3) && ~isempty(ST_STmodels3(model_idxave3).packavegrage));
    ST_model3      = []; if has_ST3,      ST_model3      = ST_STmodels3(model_idx3).packed;            end
    ST_model_avg3  = []; if has_average3, ST_model_avg3  = ST_STmodels3(model_idxave3).packavegrage;   end

    % 时间窗候选（ST 分支用过去120s；UnK 用当前历元）
    time_range    = est_time - 120 : est_time;
    time_idx_win  = find(ismember(IPPVD(:,1), time_range));
    time_idx_now  = find(IPPVD(:,1) == est_time);

    if isempty(time_idx_win)
        % ST 需要窗口；若窗口空，跳过该历元
        continue;
    end
    if isempty(time_idx_now) && ~(has_ST||has_average||has_ST2||has_average2||has_ST3||has_average3)
        % 兜底也需要当前历元 IPP；且六模型全无 → 无法估计
        continue;
    end

    % 日内 30s 索引
    j  = epoch_to_30s_index(est_time);
    jj = j - 1;  % 前一个30s

    for ig = 1:num_igps
        % Broadcast GIVD 真值
        truth_val = NaN; truth_val2 = NaN;
        m = B_map(ig);
        if ~isnan(m)
            if ~isnan(j)  && j>=1  && j<=2880,  truth_val  = Broadcast_GIVD(m, j);  end
            if ~isnan(jj) && jj>=1 && jj<=2880, truth_val2 = Broadcast_GIVD(m, jj); end
        end

        lat1 = IGP_lat(ig);
        lon1 = IGP_lon(ig);

        % 进度提示
        if ~isempty(dq) && mod(ig, 50) == 1
            task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
            send(dq, sprintf('线程 %d | epoch %d/%d | IGP %d/%d', tid, ie, num_epochs, ig, num_igps));
        end

        % ============ 统一六模型独立排序→选点→Kriging ============
        % 候选样本地理
        lat_cand = IPPVD(time_idx_win, 3);
        lon_cand = IPPVD(time_idx_win, 4);

        % 候选样本与 IGP 的大圆距离
        dlat = deg2rad(lat_cand - lat1);
        dlon = deg2rad(lon_cand - lon1);
        a    = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat_cand)).*sin(dlon/2).^2;
        c    = 2 * atan2(sqrt(a), sqrt(1-a));
        Dist_all = R_earth * c;

        selected = (Dist_all <= R_max & Dist_all > 0);
        id_all   = time_idx_win(selected);
        if numel(id_all) < N_min
            % 候选样本不足：若六模型全无 → 兜底也无意义；否则跳过该 IGP
            continue;
        end

        % stpacked 所需量
        s     = Dist_all(selected);                 % 空间距离 (km)
        tdiff = abs(IPPVD(id_all,1) - est_time);    % 时间差 (s)

        goto_UNK = ~(has_ST || has_average || has_ST2 || has_average2 || has_ST3 || has_average3);

        % best 记录器
        got_best   = false;
        best_err   = Inf;
        best_vec   = [];
        best_COND  = [NaN,NaN,NaN];
        best_Rfitt = NaN; best_R_T = NaN; best_TCM = NaN; best_RCM1 = NaN;

        % 六个候选（顺序可按偏好调整）
        cands = struct( ...
          'has',   {has_ST, has_average, has_ST2, has_average2, has_ST3, has_average3}, ...
          'model', {ST_model, ST_model_avg, ST_model2, ST_model_avg2, ST_model3, ST_model_avg3}, ...
          'label', {'ST1','avg1','ST2','avg2','ST3','avg3'} );

        for cc = 1:numel(cands)
            if ~cands(cc).has, continue; end

            % 计算该模型的 stpacked 距离→排序
            v_st       = stpacked(cands(cc).model, s, tdiff);
            [~, ord]   = sort(v_st, 'ascend');
            id_sel     = id_all(ord);
            if numel(id_sel) > N_target, id_sel = id_sel(1:N_target); end
            PF         = IPPVD(id_sel, :);
            if size(PF,1) < N_min, continue; end

            % 距离/时间量
            dlat = deg2rad(PF(:,3) - lat1);
            dlon = deg2rad(PF(:,4) - lon1);
            aa   = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(PF(:,3))).*sin(dlon/2).^2;
            ccgc = 2 * atan2(sqrt(aa), sqrt(1-aa));
            Dist = R_earth * ccgc;
            t_remaining = est_time - PF(:,1);

            % Rfitt/R_T（凑满 N_target 用 max(Dist)，否则 R_max）
            if numel(id_sel) >= N_target
                Rfitt = max(max(Dist), R_min);
            else
                Rfitt = R_max;
            end
            R_T = max(t_remaining);

            % East/North
            East2  = zeros(size(PF,1),1);
            North2 = zeros(size(PF,1),1);
            for k2 = 1:size(PF,1)
                [E2, N2, ~] = geodetic2enu(PF(k2,3), PF(k2,4), height, lat1, lon1, height, wgs84, 'degrees');
                East2(k2)  = E2;  North2(k2) = N2;
            end

            % 距离矩阵与时间差矩阵
            ns = size(PF,1);
            dist_n0 = distance('gc', PF(:,3), PF(:,4), lat1, lon1, R_earth);
            [Mg, Ng] = meshgrid(1:ns);
            dist_nn = reshape(distance('gc', ...
                PF(Mg(:),3), PF(Mg(:),4), PF(Ng(:),3), PF(Ng(:),4), R_earth), ns, ns);
            t_nn = reshape(abs(PF(Mg(:),1) - PF(Ng(:),1)), ns, ns);

            % RCM/TCM
            lati_mean  = mean(PF(:,3)); long_mean = mean(PF(:,4));
            Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
            RCM1       = Rcentroid1 / max(Rfitt, eps);
            if R_T == 0, TCM = 0; else, TCM = mean(t_remaining) / R_T; end

            % STkriging
            data_values = PF(:,5);
            [X0, Y0]    = meshgrid(0);
            [uk, COND]  = STkriging(cands(cc).model, East2, North2, t_remaining, data_values, ...
                                    X0, Y0, 0, dist_nn, dist_n0, t_nn, truth_val);
            if isempty(uk), continue; end

            % 与真值误差比较
            cur_err = abs(uk(1) - truth_val);
            if isfinite(cur_err)
                if cur_err < best_err
                    got_best   = true; best_err = cur_err;
                    best_vec   = uk;   best_COND = COND;
                    best_Rfitt = Rfitt; best_R_T = R_T; best_TCM = TCM; best_RCM1 = RCM1;
                end
            elseif ~got_best
                got_best   = true; best_err = cur_err;
                best_vec   = uk;   best_COND = COND;
                best_Rfitt = Rfitt; best_R_T = R_T; best_TCM = TCM; best_RCM1 = RCM1;
            end
        end % for cc（六候选）

        % ---------------- 写出或兜底 ----------------
        if got_best
            % 只在选出最佳候选后，计算一次 sigma_undersampled / sigma_GIVE 并写结果
            sigma_undersampled = NaN; sigma_GIVE = NaN;
            if strcmpi(type,'GIVE') && ~isempty(monothreat)
                Rfit_bin = discretize(best_Rfitt, Rfit_edges);
                RCM_bin  = discretize(best_RCM1, RCM_edges);
                TCM_bin  = discretize(best_TCM,  TCM_edges);
                if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                    sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                    if numel(best_vec) >= 6 && isfinite(best_vec(6))
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
            continue; % 本 IGP 完成

        elseif goto_UNK
            % ========== 六模型全无：进入 UnK 兜底 ==========
            if isempty(time_idx_now), continue; end
            PF_now = IPPVD(time_idx_now, :);
            if isempty(PF_now), continue; end

            % 计算所有 IPP 与 IGP 的距离
            lat2 = PF_now(:,3);  lon2 = PF_now(:,4);
            dlat = deg2rad(lat2 - lat1);
            dlon = deg2rad(lon2 - lon1);
            a = sin(dlat/2).^2 + cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(dlon/2).^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            Dist_full = R_earth * c;

            % 排序并跳过自身（最小可能是 IGP==IPP）
            [Dist_sorted, idx_sort] = sort(Dist_full, 'ascend');
            if ~isempty(Dist_sorted), Dist_sorted(1) = []; idx_sort(1) = []; end

            % 确定 R_fit
            num_IPPs_le_Rmin = sum(Dist_sorted <= R_min);
            if num_IPPs_le_Rmin >= N_target
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
            if size(PF_samples,1) < N_min, continue; end

            % East/North
            East2  = zeros(size(PF_samples,1),1);
            North2 = zeros(size(PF_samples,1),1);
            for k = 1:size(PF_samples,1)
                [E, Nn, ~] = geodetic2enu(PF_samples(k,3), PF_samples(k,4), height, ...
                    lat1, lon1, height, wgs84, 'degrees');
                East2(k)  = E;  North2(k) = Nn;
            end

            % 距离矩阵（无时间维）
            ns = size(PF_samples,1);
            dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
            [Mgrid, Ngrid] = meshgrid(1:ns);
            dist_nn = reshape(distance('gc', ...
                PF_samples(Mgrid(:),3), PF_samples(Mgrid(:),4), ...
                PF_samples(Ngrid(:),3), PF_samples(Ngrid(:),4), R_earth), ns, ns);

            % RCM（无时间维）/ R_T、TCM=0
            lati_mean  = mean(PF_samples(:,3)); long_mean = mean(PF_samples(:,4));
            Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
            RCM1 = Rcentroid1 / max(Rfitt, eps);
            R_T  = 0; TCM = 0;

            % UnK_threat
            data_values = PF_samples(:,5);
            [X0, Y0] = meshgrid(0);
            UnK_result = Unk_threat(defaultUnK, East2, North2, data_values, X0, Y0, defaultfit, dist_nn, dist_n0);
            % UnK_result: [zi, chi, chi_irreg, R_irreg, sigma_IPP, sigma_inflated]

            % sigma_undersampled / GIVE
            sigma_undersampled = NaN; sigma_GIVE = NaN;
            if ~isempty(monothreat)
                if ndims(monothreat) == 2
                    Rfit_bin = discretize(Rfitt, Rfit_edges);
                    RCM_bin  = discretize(RCM1, RCM_edges);
                    if ~isnan(Rfit_bin) && ~isnan(RCM_bin)
                        sigma_undersampled = monothreat(RCM_bin, Rfit_bin);
                    end
                else
                    Rfit_bin = discretize(Rfitt, Rfit_edges);
                    RCM_bin  = discretize(RCM1, RCM_edges);
                    TCM_bin  = 1; % 假定第一层为 TCM=0
                    if ~isnan(Rfit_bin) && ~isnan(RCM_bin) && ~isnan(TCM_bin)
                        sigma_undersampled = monothreat(TCM_bin, RCM_bin, Rfit_bin);
                    end
                end
                if isfinite(sigma_undersampled)
                    sigma_GIVE = sqrt(sigma_undersampled.^2 + UnK_result(6).^2);
                end
            end

            % 写结果
            save_row = save_row + 1;
            results(save_row, :) = [ ...
                est_time, truth_val, truth_val2, lat1, lon1, ...
                UnK_result(:).', ...
                Rfitt, RCM1, R_T, TCM, ...
                sigma_undersampled, sigma_GIVE ...
            ];
            CONDs(save_row, :) = [NaN, NaN, NaN];
            continue;

        else
            % 有模型但因样本不足/kriging失败均未得到有效结果：跳过该 IGP
            continue;
        end

    end % for ig
end % for ie

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
