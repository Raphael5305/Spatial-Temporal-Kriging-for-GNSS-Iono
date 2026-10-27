function results = STIPP_fast(ST_STmodels, ST_STmodels2, ST_STmodels3, R, N, IPPVD, samp_int, Rin_numbers, IGP, dq)
% STIPPthreat_fast_IPPbox — 六模型独立估计（无 data deprivation）
% 目标：在 IGP±2.5° 盒内，以实际观测点为 target（其第5列为“真值”），
% 在过去120 s 时间窗中，按各 ST 模型的 stpacked(s,tdiff) 排序选样、STkriging 估计，
% 六个模型分别估计，取 |估计−真值| 最小者为最佳，写出一行结果。
%
% 输出（每行 15 列）:
% [ epoch, sat, lat, lon, truth, uk_best(1:6), Rfitt, RCM1, R_T, TCM ]

if nargin < 10, dq = []; end %#ok<*NASGU>

% 参数
R_min    = R(1); % 800
R_max    = R(2); % 2100
N_min    = N(1); % 10
N_target = N(2); % 30

% 常量
wgs84   = wgs84Ellipsoid('kilometer');
height  = 350;          % km
R_earth = 6378.137;     % km
DEG2RAD = pi/180;

% 数据场
PF        = IPPVD;             % [epoch, sat, lat, lon, meas, ...]
pf_epoch  = PF(:,1);
pf_lat    = PF(:,3);
pf_lon    = mod(PF(:,4) + 180, 360) - 180;   % 规范到 [-180, 180]

% IGP（仅一行,lon优化为-180到180）
IGP_lat = IGP(1);
IGP_lon = mod(IGP(2) + 180, 360) - 180;

% 历元集合（按 samp_int 下采样）
all_epochs = unique(pf_epoch, 'stable');
est_epochs = all_epochs(1:samp_int:numel(all_epochs));
num_epochs = numel(est_epochs);

% 预分配（每个 target 一行；按经验给出保守上限）
alloc_step = max(1, Rin_numbers*9*50);  % 仅作缓冲用
results  = zeros(alloc_step, 15);
save_row = 0;

for ie = 1:num_epochs
    % 进度提示
    if ~isempty(dq) && mod(ie, 400) == 1
        t = getCurrentTask(); if ~isempty(t), tid = t.ID; else, tid = 0; end
        send(dq, sprintf('STK线程 %d | epoch %d/%d', tid, ie, num_epochs));
    end

    est_time = est_epochs(ie);

    % === 模型就位（与 NLSSTGIVx2 相同的六模型）===
    % 第一套
    model_idx    = find([ST_STmodels.epoch_time] == est_time, 1, 'first');
    model_idxave = []; if ~isempty(model_idx), model_idxave = ceil(model_idx/10)*10; end
    if ~isempty(model_idxave), model_idxave = min(2870, model_idxave); end
    has_ST       = (~isempty(model_idx)    && ~isempty(ST_STmodels(model_idx).packed));
    has_average  = (~isempty(model_idxave) && ~isempty(ST_STmodels(model_idxave).packavegrage));
    ST_model     = []; if has_ST,      ST_model     = ST_STmodels(model_idx).packed;          end
    ST_model_avg = []; if has_average, ST_model_avg = ST_STmodels(model_idxave).packavegrage; end

    % 第二套
    model_idx2    = find([ST_STmodels2.epoch_time] == est_time, 1, 'first');
    model_idxave2 = []; if ~isempty(model_idx2), model_idxave2 = ceil(model_idx2/10)*10; end
    if ~isempty(model_idxave2), model_idxave2 = min(2870, model_idxave2); end
    has_ST2       = (~isempty(model_idx2)    && ~isempty(ST_STmodels2(model_idx2).packed));
    has_average2  = (~isempty(model_idxave2) && ~isempty(ST_STmodels2(model_idxave2).packavegrage));
    ST_model2      = []; if has_ST2,      ST_model2      = ST_STmodels2(model_idx2).packed;          end
    ST_model_avg2  = []; if has_average2, ST_model_avg2  = ST_STmodels2(model_idxave2).packavegrage; end

    % 第三套
    model_idx3    = find([ST_STmodels3.epoch_time] == est_time, 1, 'first');
    model_idxave3 = []; if ~isempty(model_idx3), model_idxave3 = ceil(model_idx3/10)*10; end
    if ~isempty(model_idxave3), model_idxave3 = min(2870, model_idxave3); end
    has_ST3       = (~isempty(model_idx3)    && ~isempty(ST_STmodels3(model_idx3).packed));
    has_average3  = (~isempty(model_idxave3) && ~isempty(ST_STmodels3(model_idxave3).packavegrage));
    ST_model3      = []; if has_ST3,      ST_model3      = ST_STmodels3(model_idx3).packed;          end
    ST_model_avg3  = []; if has_average3, ST_model_avg3  = ST_STmodels3(model_idxave3).packavegrage; end

    if ~(has_ST||has_average||has_ST2||has_average2||has_ST3||has_average3)
        % 六模型全无 → 无法估计
        continue;
    end

    cands = struct( ...
        'has',   {has_ST, has_average, has_ST2, has_average2, has_ST3, has_average3}, ...
        'model', {ST_model, ST_model_avg, ST_model2, ST_model_avg2, ST_model3, ST_model_avg3}, ...
        'label', {'ST1','avg1','ST2','avg2','ST3','avg3'} );

    % 时间窗口（过去120 s）
    mask_win     = (pf_epoch >= est_time - 120) & (pf_epoch <= est_time);
    time_idx_win = find(mask_win);
    if isempty(time_idx_win), continue; end

    % 本历元且落在 IGP ±2.5° 盒内的 target 候选（以规范化经度判断）
    mask_box_exact = (pf_epoch == est_time) & ...
                     (abs(pf_lat - IGP_lat) <= 2.5) & ...
                     (abs(pf_lon - IGP_lon) <= 2.5);
    target_idx = find(mask_box_exact);
    if isempty(target_idx), continue; end

    % -------- 遍历所有 target 点 --------
    for tk = 1:numel(target_idx)
        k = target_idx(tk);

        % target 基本量（注意：经度用规范化后的 pf_lon，但 target_point 保留原 PF 第1~5列）
        lat1 = pf_lat(k);
        lon1 = pf_lon(k);
        truth_val = PF(k,5);     % 第5列作为“真值”

        % 候选样本地理与大圆距离（相对 target）
        lat_cand = pf_lat(time_idx_win);
        lon_cand = pf_lon(time_idx_win);

        dlat = (lat_cand - lat1)*DEG2RAD;
        dlon = (lon_cand - lon1)*DEG2RAD;
        a    = sin(0.5*dlat).^2 + cos(lat1*DEG2RAD).*cos(lat_cand*DEG2RAD).*sin(0.5*dlon).^2;
        Dist_all = R_earth .* (2*atan2(sqrt(a), sqrt(1-a)));

        % 去除 data deprivation：**不剔除目标点**，R_max 以内全纳入
        selected = (Dist_all <= R_max&Dist_all > 0);
        id_all   = time_idx_win(selected);

        if numel(id_all) < N_min
            % 候选样本不足，跳过该 target
            continue;
        end

        % stpacked 所需量
        s     = Dist_all(selected);                   % 空间距离 (km)
        tdiff = abs(pf_epoch(id_all) - est_time);     % 时间差 (s)

        % 记录器
        got_best  = false;
        best_err  = Inf;
        best_vec  = [];
        best_Rfitt= NaN; best_R_T = NaN; best_TCM = NaN; best_RCM1 = NaN;

        % 六个候选（按可用性循环）
        for cc = 1:numel(cands)
            if ~cands(cc).has, continue; end

            % 该模型的时空“距离”→排序
            v_st     = stpacked(cands(cc).model, s, tdiff);
            [~, ord] = sort(v_st, 'ascend');
            id_sel   = id_all(ord);
            if numel(id_sel) > N_target, id_sel = id_sel(1:N_target); end
            PFsel    = PF(id_sel, :);
            if size(PFsel,1) < N_min, continue; end

            % 量化 Rfitt / R_T
            dlat2 = (PFsel(:,3) - lat1)*DEG2RAD;
            dlon2 = (PFsel(:,4) - lon1)*DEG2RAD;
            aa    = sin(0.5*dlat2).^2 + cos(lat1*DEG2RAD).*cos(PFsel(:,3)*DEG2RAD).*sin(0.5*dlon2).^2;
            Dist  = R_earth .* (2*atan2(sqrt(aa), sqrt(1-aa)));
            t_remaining = est_time - PFsel(:,1);

            if numel(id_sel) >= N_target
                Rfitt = max(max(Dist), R_min);
            else
                Rfitt = R_max;
            end
            R_T = max(t_remaining);

            % East/North（相对 target 位置）
            East2  = zeros(size(PFsel,1),1);
            North2 = zeros(size(PFsel,1),1);
            for k2 = 1:size(PFsel,1)
                [E2, N2, ~] = geodetic2enu(PFsel(k2,3), PFsel(k2,4), height, lat1, lon1, height, wgs84, 'degrees');
                East2(k2)  = E2;  North2(k2) = N2;
            end

            % 构造距离/时间差矩阵
            ns = size(PFsel,1);
            dist_n0 = distance('gc', PFsel(:,3), PFsel(:,4), lat1, lon1, R_earth);
            [Mg, Ng] = meshgrid(1:ns);
            dist_nn = reshape(distance('gc', ...
                PFsel(Mg(:),3), PFsel(Mg(:),4), PFsel(Ng(:),3), PFsel(Ng(:),4), R_earth), ns, ns);
            t_nn = reshape(abs(PFsel(Mg(:),1) - PFsel(Ng(:),1)), ns, ns);

            % RCM/TCM
            lati_mean  = mean(PFsel(:,3)); long_mean = mean(PFsel(:,4));
            Rcentroid1 = distance('gc', lati_mean, long_mean, lat1, lon1, R_earth);
            RCM1       = Rcentroid1 / max(Rfitt, eps);
            if R_T == 0, TCM = 0; else, TCM = mean(t_remaining) / R_T; end

            % STkriging 估计
            data_values = PFsel(:,5);
            [X0, Y0]    = meshgrid(0);
            [uk, ~]     = STkriging(cands(cc).model, East2, North2, t_remaining, data_values, ...
                                    X0, Y0, 0, dist_nn, dist_n0, t_nn, truth_val);
            if isempty(uk), continue; end

            cur_err = abs(uk(1) - truth_val);
            if isfinite(cur_err) && cur_err < best_err
                got_best  = true; best_err = cur_err;
                best_vec  = uk;
                best_Rfitt= Rfitt; best_R_T = R_T; best_TCM = TCM; best_RCM1 = RCM1;
            elseif ~got_best && ~isfinite(best_err)
                % 作为首次可用
                got_best  = true; best_err = cur_err;
                best_vec  = uk;
                best_Rfitt= Rfitt; best_R_T = R_T; best_TCM = TCM; best_RCM1 = RCM1;
            end
        end % for cc

        if ~got_best
            % 六模型都失败：跳过该 target
            continue;
        end

        % 写结果（如需 16 列，可在末尾加一个 NaN 占位）
        target_point = PF(k, 1:5);
        row = [ target_point, best_vec(:).', best_Rfitt, best_RCM1, best_R_T, best_TCM ];

        if save_row >= size(results,1)
            results = [results; zeros(alloc_step, size(results,2))]; %#ok<AGROW>
        end
        save_row = save_row + 1;
        results(save_row, :) = row;

    end % target loop
end % epoch loop

% 裁剪
results = results(1:save_row, :);
end
