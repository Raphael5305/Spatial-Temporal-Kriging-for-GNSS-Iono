function results = STIPPthreat_fast_IPPbox(ST_STmodels, ST_STmodels2,ST_STmodels3,R, N, IPPVD, samp_int, Rin_numbers, IGP, dq)
% 基于 IPPVD 的局部框选 + 多 target 点逐个剥离
if nargin <10, dq = []; end

R_min    = R(1); % 800
R_max    = R(2); % 2100
N_min    = N(1); % 10
N_target = N(2); % 30

% 常量
wgs84   = wgs84Ellipsoid('kilometer'); %#ok<NASGU>
height  = 350;          % km
R_earth = 6378.137;     % km
DEG2RAD = pi/180;

PF        = IPPVD; % [epoch, sat, lat, lon, meas, ...]
pf_epoch  = PF(:,1);
pf_lat    = PF(:,3);
pf_lon    = mod(PF(:,4) + 180, 360) - 180;   % 仅用于距离与筛选的经度规范

% IGP（仅一行）
IGP_lat = IGP(1);
IGP_lon = mod(IGP(2) + 180, 360) - 180;

% 历元集合（按 samp_int 下采样）
all_epochs = unique(pf_epoch, 'stable');
est_epochs = all_epochs(1:samp_int:numel(all_epochs));
num_epochs = numel(est_epochs);

% 预分配并按需扩容
alloc_step = max(1, Rin_numbers*9*5*50*2880);
results  = zeros(alloc_step, 16);
save_row = 0;

for ie = 1:num_epochs
    % 进度提示
    if ~isempty(dq) && mod(ie, 400) == 1
        t = getCurrentTask(); if ~isempty(t), tid = t.ID; else, tid = 0; end
        send(dq, sprintf('STK线程 %d | epoch %d/%d', tid, ie, num_epochs));
    end

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

    if  ~(has_ST||has_average||has_ST2||has_average2||has_ST3||has_average3)
        % 六模型全无 → 无法估计
        continue;
    end

    cands = struct( ...
        'has',   {has_ST, has_average, has_ST2, has_average2, has_ST3, has_average3}, ...
        'model', {ST_model, ST_model_avg, ST_model2, ST_model_avg2, ST_model3, ST_model_avg3});

% 时间窗口
    mask_win      = (pf_epoch >= est_time-120) & (pf_epoch <= est_time);
    time_idx_win  = find(mask_win);

    % 本历元且落在 IGP ±2.5° 盒内的 target 候选（按规范化经度判断）
    mask_box_exact = (pf_epoch == est_time) & ...
                     (abs(pf_lat - IGP_lat) <= 2.5) & ...
                     (abs(pf_lon - IGP_lon) <= 2.5);
    target_idx = find(mask_box_exact);
    if isempty(target_idx)
        continue;
    end



    % 遍历所有 target 点
    for tk = 1:numel(target_idx)
        k = target_idx(tk);

        % 注意：距离用规范化后的经度；target_point 则直接取原 PF 的前5列
        lat1 = pf_lat(k);
        lon1 = pf_lon(k);
        cos_lat1 = cos(lat1*DEG2RAD);

        % == 你要求的改动：target_point 直接存该点的前5列（第5列为测量值=真值） ==
        target_point = PF(k, 1:5);

        % 在时间窗内挑样本（剔除与 target 重合点）
        lat_cand = pf_lat(time_idx_win);
        lon_cand = pf_lon(time_idx_win);

        % 大圆距离
        dlat = (lat_cand - lat1)*DEG2RAD;
        dlon = (lon_cand - lon1)*DEG2RAD;
        a = sin(0.5*dlat).^2 + cos_lat1 .* cos(lat_cand*DEG2RAD) .* sin(0.5*dlon).^2;
        Dist_all = R_earth .* (2*atan2(sqrt(a), sqrt(1-a)));

        selected = (Dist_all <= R_max & Dist_all > 0);
        idsearch_initial = time_idx_win(selected);
        if numel(idsearch_initial) < N_min
            continue;
        end

        % ST 变差函数排序（升序）
        s     = Dist_all(selected);
        tdiff = abs(pf_epoch(idsearch_initial) - est_time);

        mask_model = [];
        if  ~isempty(cands(5).model)
            mask_model = cands(5).model;   % ST3
        else
            % 回退：找第一个可用的
            for ii = 1:numel(cands)
                if cands(ii).has && ~isempty(cands(ii).model)
                    mask_model = cands(ii).model;
                    break;
                end
            end
        end
        vvals = stpacked(mask_model, s, tdiff);
        [sorted_variogram, ord] = sort(vvals, 'ascend');
        idsearch_initial = idsearch_initial(ord);
        variogram_values = sorted_variogram;

        max_variogram = max(variogram_values);
        min_variogram = min(variogram_values);
        PF_samples    = PF(idsearch_initial, :);

        if numel(idsearch_initial) >= N_target
            variogram_des = variogram_values(N_target);
        else
            variogram_des = max_variogram;
        end

        % -------- Mode 1: Fixed width --------
        R_in_values = linspace(min_variogram-1e-10, variogram_des+1e-10, Rin_numbers+1);
        width_all   =[0, (0.1:0.2:0.9) .* (max_variogram - min_variogram)];
        for r = 1:length(R_in_values)-1
            R_in_left = R_in_values(r);
            for Rw = width_all
                R_in  = R_in_left;
                R_out = R_in_left + Rw;
                result_once =STdepriveAndEstimate2_multimodels(R_min, R_max, N, R_in, R_out, ...
                                PF_samples, target_point, cands, wgs84, height, R_earth, 'fixed width');
                if isnan(result_once), continue; end
                if save_row >= size(results,1)
                    results = [results; zeros(alloc_step, size(results,2))]; %#ok<AGROW>
                end
                save_row = save_row + 1;
                results(save_row, :) = [target_point result_once];
            end
        end

        % -------- Mode 3: Temporal delete --------
        R_in_values2 = 0:30:120;
        Rw_values2   = (0.1:0.2:0.9) * 120;
        for r = 1:length(R_in_values2)
            R_in = R_in_values2(r);
            for w = 1:numel(Rw_values2)
                R_out = R_in + Rw_values2(w);
                result_once =STdepriveAndEstimate2_multimodels(R_min, R_max, N, R_in, R_out, ...
                                PF_samples, target_point, cands, wgs84, height, R_earth, 'temporal delete');
                if isnan(result_once), continue; end
                if save_row >= size(results,1)
                    results = [results; zeros(alloc_step, size(results,2))]; %#ok<AGROW>
                end
                save_row = save_row + 1;
                results(save_row, :) = [target_point result_once];
            end
        end

    end % target loop
end % epoch loop

% 裁剪
results = results(1:save_row, :);
end
