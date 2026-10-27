function save_IGP_variogram(PF, igp_lat, igp_lon, params, outdir)
% 保存“单个IGP”的每一次历元经验变异函数与拟合模型到一个文件：
%   文件名：IGP_<latNS>_<lonEW>_variogram.mat   例：IGP_30N_100W_variogram.mat
%
% 输入:
%   PF      : [epoch, sat_id, lat, lon, y, ...]  其中 lon ∈ [-180,180]（保留原样）
%   igp_lat : IGP中心纬度（例如 30）
%   igp_lon : IGP中心经度（例如 -100 表示 100°W）
%   params  : 结构体，示例字段：
%             .sth_km = 400;       % 空间邻近阈值 (km)
%             .win_deg = 5;        % IGP矩形窗口半宽 (deg)
%             .time_back = 300;    % 时间滑窗长度 (s)，10个历元
%             .time_step = 30;     % 历元步长 (s)
%             .use_single = true;  % 存盘前是否压成 single
%             .code_version = 'v1.0';
%   outdir  : 输出目录
%
% 输出:
%   IGP, epochs, params_saved, ST_S, ST_T, ST_ST, ST_STmodels（见下）

%% Example：
% params = struct('sth_km',400,'win_deg',5,'time_back',150,'time_step',30,...
%                 'use_single',true,'code_version','v1.2');
% outdir = 'E:\CONUS_IGP_variograms';
% [LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
% IGP_list=[LatC(:), LonC(:)];
% for i = 1:size(IGP_list,1)
%     save_IGP_variogram(PF, IGP_list(i,1), IGP_list(i,2), params, outdir);
% end
%% -------------------------------------------------------------

    if nargin < 5, outdir = pwd; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % —— 参数默认值 ——
    p = struct('sth_km',2100,'win_deg',15,'time_back',300,'time_step',30,...
               'use_single',true,'code_version','v1.0');

    if exist('params','var') && ~isempty(params)
        fns = fieldnames(params);
        for i=1:numel(fns), p.(fns{i}) = params.(fns{i}); end
    end

    % —— 历元集合 ——
    all_epochs = unique(PF(:,1));

    % —— 预分配 cell/struct ——
    % ST_S  = cell(numel(all_epochs),1);
    % ST_T  = cell(numel(all_epochs),1);
    % ST_ST = cell(numel(all_epochs),1);
    ST_STmodels = struct('epoch_time',cell(numel(all_epochs),1), ...
                         'model',cell(numel(all_epochs),1), ...
                         'Smodel',cell(numel(all_epochs),1), ...
                         'Tmodel',cell(numel(all_epochs),1));
    kept = false(numel(all_epochs),1);

    % —— 主循环：逐历元在 IGP 局部统计 + 拟合 ——
    for k = 1:numel(all_epochs)
        t = all_epochs(k); %当前历元：s

        in_window = PF(:,1) >= (t - p.time_back) & PF(:,1) <= t;
        if ~any(in_window), continue; end
        PF_t = PF(in_window, :);

        % 进一步IGP 矩形窗口（经度使用 -180..180，允许跨 -180/180）
        roi = PF_t(:,3) >= (igp_lat - p.win_deg) & PF_t(:,3) <= (igp_lat + p.win_deg) & ...
              in_lon_window_180(PF_t(:,4), igp_lon, p.win_deg);
        if nnz(roi) < 5, continue; end

        % 构造 x_sub / y_sub（y 假定在第5列，按你的数据修改）
        x_sub = [PF_t(roi,3), PF_t(roi,4), PF_t(roi,1)];  % [lat, lon(-180..180), epoch]
        y_sub =  PF_t(roi,7); %residuals

        ST_S_loc  = variogram_spatial([x_sub(:,1:2), x_sub(:,3)], y_sub, 'sth', p.sth_km);
        ST_T_loc  = variogram_temporal2([x_sub(:,1:2), x_sub(:,3), y_sub], 'sth', 15,'n',p.time_back/p.time_step);
        ST_ST_loc = variogram_ST2([x_sub(:,1:2), x_sub(:,3)], y_sub, 'sth', p.sth_km, 'n', p.time_back/p.time_step);

        if isempty(ST_S_loc) && isempty(ST_T_loc) && isempty(ST_ST_loc), continue; end 
        %若均为空不拟合ST函数，在真实处理时无ST考虑使用该处的S函数退化为universal kriging或者CONUS的unk

        % 拟合（与你 STthreat 的 Model 2 一致）
        try
            [ST2, Smodel, Tmodel] = fit_ST2_from_loc(ST_S_loc, ST_T_loc, ST_ST_loc);
        catch ME
            warning('Epoch %d: model fitting failed: %s', t, ME.message);
            ST2 = []; Smodel = []; Tmodel = [];
        end

        % if p.use_single
        %     ST_S{k}  = to_single(ST_S_loc);
        %     ST_T{k}  = to_single(ST_T_loc);
        %     ST_ST{k} = to_single(ST_ST_loc);
        % else
        %     ST_S{k}  = ST_S_loc;
        %     ST_T{k}  = ST_T_loc;
        %     ST_ST{k} = ST_ST_loc;
        % end

        ST_STmodels(k).epoch_time = t;
        ST_STmodels(k).model  = ST2;
        % ST_STmodels(k).Smodel = Smodel;
        % ST_STmodels(k).Tmodel = Tmodel;

        kept(k) = true;
    end

    % —— 只保留有效历元 ——
    % ST_S  = ST_S(kept);
    % ST_T  = ST_T(kept);
    % ST_ST = ST_ST(kept);
    ST_STmodels = ST_STmodels(kept);
    epochs = all_epochs(kept);

    % —— 落盘 ——
    IGP.lat = igp_lat;
    IGP.lon = igp_lon;
    params_saved = p;

    outfile = fullfile(outdir, sprintf('IGP_%s_%s_variogram.mat', tagLatNS(igp_lat), tagLonEW(igp_lon)));
    save(outfile, 'IGP', 'epochs', 'params_saved', 'ST_STmodels', '-v7.3');
    % save(outfile, 'IGP', 'epochs', 'params_saved', 'ST_S', 'ST_T', 'ST_ST', 'ST_STmodels', '-v7.3');
    fprintf('[OK] Saved %s (epochs=%d)\n', outfile, numel(epochs));
end

% ---------- 辅助函数们 ----------

function tf = in_lon_window_180(lon, lon0, halfwin_deg)
% lon, lon0 ∈ [-180,180]，判断 lon 是否落在以 lon0 为中心、半宽 halfwin_deg 的经度窗口
% 采用最短角距：d = wrapTo180(lon - lon0); |d| <= halfwin
    d = lon - lon0;
    d = mod(d + 180, 360) - 180;  % wrap 到 [-180,180]
    tf = abs(d) <= halfwin_deg;
end

function Xs = to_single(X)
    if istable(X)
        Xs = single(table2array(X));
    else
        Xs = single(X);
    end
end

function [ST2, Smodel, Tmodel] = fit_ST2_from_loc(ST_S, ST_T, ST_ST)
    if istable(ST_T), ST_T = table2array(ST_T); end
    if istable(ST_S), ST_S = table2array(ST_S); end
    if istable(ST_ST), ST_ST = table2array(ST_ST); end

    % 空间
    if ~isempty(ST_S)
        s_lags_s = ST_S(:,1); exp_values_s = ST_S(:,2);
        [~, s_centers, s_grouped_values] = bin_grouping(s_lags_s, exp_values_s, 50);
        Smodel = PSO_maternfit(s_centers, s_grouped_values, 0.01, 0.8);
    else
        Smodel = [];
    end

    % 时间,考虑有4个时间差对以上的情况，否则难以解算PSO模型
    if ~isempty(ST_T)
        t_lags_t = ST_T(:,1); exp_values_t = ST_T(:,2);
        [t_unique, ~, t_group_idx] = unique(t_lags_t);
        t_grouped_values = accumarray(t_group_idx, exp_values_t, [], @mean);
        Tmodel = PSO_maternfit(t_unique, t_grouped_values, 0.01, 0.8);
    else
        Tmodel = [];
    end

    % 时空 + k1,k2,k3
    if isempty(ST_ST) || isempty(Smodel) || isempty(Tmodel)
        ST2 = []; return;
    end

    s_lags = ST_ST(:,1);
    t_lags = ST_ST(:,2);
    exp_values = ST_ST(:,3);

    num_s_bins = 50;
    [s_centers2, t_unique2, grouped_grid] = bin_grouping_2d(s_lags, t_lags, exp_values, num_s_bins);
    %grouped_grid为第i个时间，第j个空间的二维组，这与下面的mesh应该保持相反
    [s_grid, t_grid] = meshgrid(s_centers2, t_unique2);

    ys0 = Smodel.func(s_grid(:));
    y0t = Tmodel.func(t_grid(:));
    Y  = grouped_grid(:);
    valid = isfinite(ys0) & isfinite(y0t) & isfinite(Y);
    ys0 = ys0(valid); y0t = y0t(valid); Y = Y(valid); %一列一列展开，因此是一组时间一组时间排列

    %参数回归，最小二乘解析解：Y=ys0+y0t-α*yS.*yT + ε
    x=ys0.*y0t;
    r=Y-ys0-y0t;
    mask = (x > 0) & isfinite(r); 
    if ~any(mask)
        % 若没有可用信息点，回退 α=0（无交互）
        alpha_hat = 0;
    else
        % 加权（若你有每格 pair 数，可放到 w 里；这里默认全 1）
        w = ones(sum(mask),1);
        xm = x(mask); rm = r(mask);

        % 闭式解：alpha = sum(w*x*(yS+yT - Y)) / sum(w*x^2) = sum(-w*x*r)/sum(w*x^2)
        num = sum( w .* xm .* (-rm) );
        den = sum( w .* xm.^2 );
        if den <= eps
            alpha_hat = 0;
        else
            alpha_hat = num / den;
        end
    end

    sS=Smodel.sill;%all sill
    sT=Tmodel.sill;
    alpha_hat = max(0, min(alpha_hat, min(1/max(sS,eps), 1/max(sT,eps))));

    %组装ST2
    ST2.Smodel = Smodel;
    ST2.Tmodel = Tmodel;
    ST2.alpha  = alpha_hat;
    ST2.func   = @(s,t) ST2.Smodel.func(s)+ST2.Tmodel.func(t)-ST2.alpha.*ST2.Smodel.func(s).*ST2.Tmodel.func(t);
    ST2.nugget = Smodel.c0 + Tmodel.c0 -ST2.alpha*Smodel.c0.*Tmodel.c0;
    ST2.sill   = sS + sT - ST2.alpha * sS * sT;
end

function tag = tagLatNS(lat)
    if lat>=0, tag = sprintf('%dN', round(abs(lat)));
    else,      tag = sprintf('%dS', round(abs(lat)));
    end
end

function tag = tagLonEW(lon)
    if lon>=0, tag = sprintf('%dE', round(abs(lon)));
    else,      tag = sprintf('%dW', round(abs(lon)));
    end
end
