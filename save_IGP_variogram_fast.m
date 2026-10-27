function save_IGP_variogram_fast(PF, igp_lat, igp_lon, params, outdir)
% 保存“单个IGP”的每一次历元经验变异函数（S/T/ST）拟合结果到 MAT：
%   文件名：IGP_<latNS>_<lonEW>_variogram.mat   例：IGP_30N_100W_variogram.mat
%
% 输入:
%   PF      : [epoch, sat_id, lat, lon, ..., y] （其中 y=残差/目标量，lon∈[-180,180]）
%   igp_lat : IGP中心纬度（例如 30）
%   igp_lon : IGP中心经度（例如 -100 表示 100°W）
%   params  : 结构体，示例字段：
%             .sth_km   (默认2100)   空间邻近阈值 (km)
%             .win_deg  (默认15)     IGP矩形窗口半宽 (deg)
%             .time_back(默认300 s)  时间滑窗长度
%             .time_step(默认30 s)   历元步长
%             .use_single(默认true)  是否在存盘前压成 single（对大型格点/LUT有用）
%             .code_version(默认'v1.0')
%             .precompute_LUT (默认true) 是否预计算S/T的一维查表
%             .s_grid (可选) LUT空间距离网格，默认 0:10:sth_km
%             .t_grid (可选) LUT时间差网格，默认 0:time_step:time_back
%   outdir  : 输出目录
%
% 输出:
%   MAT内对象：IGP, epochs, params_saved, ST_STmodels
%   其中 ST_STmodels(i).packed 只含 “参数 + 可选LUT”，无函数句柄，后续评估时不跳外部函数。
    if nargin < 5, outdir = pwd; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % ------- 参数默认值 -------
    p = struct('sth_km',2100, 'win_deg',15, 'time_back',300, 'time_step',30, ...
               'use_single',true, 'code_version','v1.0', 'precompute_LUT', true);
    if exist('params','var') && ~isempty(params)
        fns = fieldnames(params);
        for i=1:numel(fns), p.(fns{i}) = params.(fns{i}); end
    end
    if ~isfield(p,'s_grid') || isempty(p.s_grid),p.s_grid = [0, 0.01, 10:10:p.sth_km];end
    if ~isfield(p,'t_grid') || isempty(p.t_grid),  p.t_grid = [0, 0.01, max(1,p.time_step):max(1,p.time_step):p.time_back]; end

    % ------- 历元集合（固定网格） -------
    full_epochs = (150:p.time_step:86400).';   % 列向量
    N = numel(full_epochs);

    % ------- 预分配（每个历元占位） -------
    ST_STmodels = struct('epoch_time',cell(N,1), ...
                         'packed',      cell(N,1), ...
                         'packavegrage',cell(N,1));
    valid_ep = false(N,1);

    % ====== 滚动块缓存（最多10条） ======
    S_blk = nan(10,4);           % [rho, nu, sigma2(partialsill), c0]
    T_blk = nan(10,4);
    blk_valid_mask = false(10,1);
    last_STST_loc  = [];         % 块内“最后一个有效历元”的 ST_ST_loc（瘦身+single 化）

    % ================= 主循环 =================
    for k = 1:N
        t = full_epochs(k);
        ST_STmodels(k).epoch_time = t;

        % 块内位置 1..10；新块开始 → 清空缓存
        pos = mod(k-1, 10) + 1;
        if pos == 1
            S_blk(:) = nan; T_blk(:) = nan;
            blk_valid_mask(:) = false;
            last_STST_loc = [];
        end

        % --------- 本历元尝试拟合（失败就保持缓存不变） ---------
        % 1) 时间窗口
        in_window = (PF(:,1) >= (t - p.time_back)) & (PF(:,1) <= t);

        if any(in_window)
            PF_t = PF(in_window, :);

            % 2) IGP 经纬窗口（lon wrap 到 [-180,180]）
            lat_ok = (PF_t(:,3) >= (igp_lat - p.win_deg)) & (PF_t(:,3) <= (igp_lat + p.win_deg));
            dlon  = PF_t(:,4) - igp_lon;
            dlon  = mod(dlon + 180, 360) - 180;
            lon_ok = abs(dlon) <= p.win_deg;
            roi = lat_ok & lon_ok;

            if nnz(roi) >= 5
                % 3) 子样本与目标量
                x_sub = [PF_t(roi,3), PF_t(roi,4), PF_t(roi,1)];  % [lat, lon, epoch]
                y_sub =  PF_t(roi,7);                             % 目标量（残差等）

                % 4) 经验变异函数
                ST_S_loc  = variogram_spatial_fast([x_sub(:,1:2), x_sub(:,3)], y_sub, 'sth', p.sth_km);
                ST_T_loc  = variogram_temporal2_fast([x_sub(:,1:2), x_sub(:,3), y_sub], 'sth', 15, 'n', p.time_back/p.time_step);
                ST_ST_loc = variogram_ST2_fast([x_sub(:,1:2), x_sub(:,3)], y_sub, 'sth', p.sth_km, 'n', p.time_back/p.time_step);

                if ~(isempty(ST_S_loc) && isempty(ST_T_loc) && isempty(ST_ST_loc))
                    % 5) 单历元拟合
                    ST2_ep = []; Smodel_ep = []; Tmodel_ep = [];
                    [ST2_ep, Smodel_ep, Tmodel_ep] = fit_ST2_from_loc(ST_S_loc, ST_T_loc, ST_ST_loc);

                    if ~isempty(ST2_ep)
                        % 标记有效
                        valid_ep(k) = true;

                        % 打包单历元
                        pk = pack_ST2(ST2_ep, p);
                        if p.use_single, pk = cast_packed_to_single(pk); end
                        ST_STmodels(k).packed = pk;

                        % 写入块缓存
                        S_blk(pos,:) = [Smodel_ep.rho, Smodel_ep.nu, Smodel_ep.sigma2, Smodel_ep.c0];
                        T_blk(pos,:) = [Tmodel_ep.rho, Tmodel_ep.nu, Tmodel_ep.sigma2, Tmodel_ep.c0];
                        blk_valid_mask(pos) = true;

                        % 更新“块内最后一个有效”的 ST_ST_loc（瘦身 + single）
                        if ~isempty(ST_ST_loc)
                            last_STST_loc = ST_ST_loc;   % 转成 single，减少内存
                        end
                    end
                end
            end
        end

        % --------- 仅此一处：块末统一处理 ---------
        if pos == 10
            if any(blk_valid_mask) && ~isempty(last_STST_loc)
                S_mean = nanmean(S_blk(blk_valid_mask,:), 1);
                T_mean = nanmean(T_blk(blk_valid_mask,:), 1);

                if all(isfinite(S_mean)) && all(isfinite(T_mean))
                    Savg = struct('rho',S_mean(1),'nu',S_mean(2),'sigma2',S_mean(3),'c0',S_mean(4));
                    Savg.sill = Savg.sigma2 + Savg.c0;
                    Savg.func = @(h) matern_eval_inline(h, struct('rho',Savg.rho,'nu',Savg.nu,'partialsill',Savg.sigma2,'nugget',Savg.c0));

                    Tavg = struct('rho',T_mean(1),'nu',T_mean(2),'sigma2',T_mean(3),'c0',T_mean(4));
                    Tavg.sill = Tavg.sigma2 + Tavg.c0;
                    Tavg.func = @(tt) matern_eval_inline(tt, struct('rho',Tavg.rho,'nu',Tavg.nu,'partialsill',Tavg.sigma2,'nugget',Tavg.c0));
                    ST2_blk = fit_ST2_alpha_from_ST_only(last_STST_loc, Savg, Tavg);
                        pk_blk  = pack_ST2(ST2_blk, p);
                        if p.use_single, pk_blk = cast_packed_to_single(pk_blk); end
                        if ~isfield(pk_blk,'meta') || isempty(pk_blk.meta), pk_blk.meta = struct; end
                        pk_blk.meta.avg_block_start_epoch = full_epochs(k-9);
                        pk_blk.meta.avg_block_end_epoch   = full_epochs(k);
                        pk_blk.meta.avg_block_size        = nnz(blk_valid_mask);
                        ST_STmodels(k).packavegrage = pk_blk;
                else
                    ST_STmodels(k).packavegrage = [];
                end
            else
                ST_STmodels(k).packavegrage = [];
            end
        end
    end

    % ------- 落盘 -------
    params_saved = p;
    epochs_all   = full_epochs;
    epochs_valid = full_epochs(valid_ep);

    outfile = fullfile(outdir, sprintf('IGP_%s_%s_variogram.mat', tagLatNS(igp_lat), tagLonEW(igp_lon)));
    try
        save(outfile,'params_saved', 'ST_STmodels','epochs_all','epochs_valid','-v7');
    catch
        save(outfile,'params_saved', 'ST_STmodels','epochs_all','epochs_valid','-v7.3');
    end
    fprintf('[OK] Saved %s (epochs_all=%d, valid=%d)\n', outfile, numel(epochs_all), numel(epochs_valid));
end
%%--------------辅助函数--------------------

function packed = pack_ST2(ST2, p)
% 从 ST2（含 Smodel, Tmodel, alpha）提取必要数值参数，并可选预计算 LUT
S = ST2.Smodel;  T = ST2.Tmodel;

% 兼容不同字段命名
S_rho = S.rho;
S_nu  = S.nu;
S_ps  = getfield_fallback(S, {'sigma2','partialsill'}); % partial sill
S_c0  = getfield_fallback(S, {'c0','nugget'});          % nugget
S_sill = getfield_fallback(S, {'sill'}, S_ps + S_c0);   % 若无 sill，则 ps+nugget

T_rho = T.rho;
T_nu  = T.nu;
T_ps  = getfield_fallback(T, {'sigma2','partialsill'});
T_c0  = getfield_fallback(T, {'c0','nugget'});
T_sill = getfield_fallback(T, {'sill'}, T_ps + T_c0);

packed.alpha  = ST2.alpha;
packed.nugget = S_c0 + T_c0 - ST2.alpha * S_c0 * T_c0;
packed.sill   = S_sill + T_sill - ST2.alpha * S_sill * T_sill;
packed.partialsill = packed.sill - packed.nugget;

% 存 S/T 参数（供重建/评估）
packed.Sparam = struct('rho',S_rho,'nu',S_nu,'partialsill',S_ps,'nugget',S_c0,'sill',S_sill);
packed.Tparam = struct('rho',T_rho,'nu',T_nu,'partialsill',T_ps,'nugget',T_c0,'sill',T_sill);

% 可选：预计算一维 LUT
if isfield(p,'precompute_LUT') && p.precompute_LUT
    packed.s_grid = p.s_grid(:).';
    packed.t_grid = p.t_grid(:).';
    packed.S_on_s = matern_eval_inline(packed.s_grid, packed.Sparam);
    packed.T_on_t = matern_eval_inline(packed.t_grid, packed.Tparam);
else
    packed.s_grid = [];
    packed.t_grid = [];
    packed.S_on_s = [];
    packed.T_on_t = [];
end

% 元信息
packed.meta = struct('code_version','packST2_v1', ...
    'win_deg', p.win_deg, ...
    'sth_km', p.sth_km, ...
    'time_back', p.time_back, ...
    'time_step', p.time_step);
end


function pk = cast_packed_to_single(pk)
% 将较大的数组压成 single，标量参数仍保留 double 的也可
f = @(x) cast(x,'single');
if ~isempty(pk.s_grid), pk.s_grid = f(pk.s_grid); end
if ~isempty(pk.t_grid), pk.t_grid = f(pk.t_grid); end
if ~isempty(pk.S_on_s), pk.S_on_s = f(pk.S_on_s); end
if ~isempty(pk.T_on_t), pk.T_on_t = f(pk.T_on_t); end
% 参数标量保持 double 以免累积误差放大（你也可以改为 single）
end

function val = getfield_fallback(S, names, default_val)
% 依次尝试多个字段名，取第一个存在的；可给默认值
for i=1:numel(names)
    if isfield(S, names{i}) && ~isempty(S.(names{i}))
        val = S.(names{i});
        return;
    end
end
if nargin >= 3
    val = default_val;
else
    error('Missing expected fields: %s', strjoin(names, '/'));
end
end

% ---------- 读取后的重建：快速 evaluator（优先 LUT） ----------
function ST_eval = rebuild_ST_eval_fast(packed)
% 基于打包内的 LUT 进行查表+线性插值；若 LUT 为空，自动降级为通用内联求值
if ~isempty(packed.s_grid) && ~isempty(packed.t_grid) ...
        && ~isempty(packed.S_on_s) && ~isempty(packed.T_on_t)
    a  = packed.alpha;
    sg = packed.s_grid; Ss = packed.S_on_s;
    tg = packed.t_grid; Tt = packed.T_on_t;
    ST_eval = @(s,t) st_eval_lut(s,t,sg,Ss,tg,Tt,a);
else
    ST_eval = rebuild_ST_eval(packed); % 无 LUT → 内联求值
end
end

function y = st_eval_lut(s,t,sg,Ss,tg,Tt,a)
% 线性插值（支持标量/向量/矩阵），边界外做线性外推
ys = interp1(sg, Ss, s, 'linear', 'extrap');
yt = interp1(tg, Tt, t, 'linear', 'extrap');
y  = ys + yt - a .* ys .* yt;       % 一次评估，不重复调用
end

function ST_eval = rebuild_ST_eval(packed)
% 无 LUT 的通用内联求值（仍不依赖 PSO_matern*）
Sparam = packed.Sparam; Tparam = packed.Tparam; a = packed.alpha;
ST_eval = @(s,t) st_eval_core(s,t,Sparam,Tparam,a);
end

function y = st_eval_core(s,t,Sparam,Tparam,a)
ys = matern_eval_inline(s, Sparam);
yt = matern_eval_inline(t, Tparam);
y  = ys + yt - a .* ys .* yt;
end

%%内联matern model function
function g = matern_eval_inline(d, P)
% 与你的 matern_variogram(d,c0,sigma2,rho,nu) 一致的半变异函数实现
% P 字段：P.rho(尺度)  P.nu(平滑度)  P.partialsill(=sigma2)  P.nugget(=c0)  P.sill(可忽略)
% --- 防护 ---
d   = max(d, 0);
rho = max(P.rho, eps);
r   = d ./ rho;                 % 无量纲距离
nu  = P.nu;
s2  = P.partialsill;            % sigma2
c0  = P.nugget;                 % c0

g = zeros(size(d), 'like', d);  % 先全置 0（确保 d=0 → γ=0）
mask = (d > 0);
if ~any(mask), return; end

rr = r(mask);

% 计算相关函数 rho_nu(rr)
if abs(nu-0.5) < 1e-12
    rho_mt = exp(-rr);
elseif abs(nu-1.5) < 1e-12
    c = sqrt(3);
    rho_mt = (1 + c*rr) .* exp(-c*rr);
elseif abs(nu-2.5) < 1e-12
    c = sqrt(5);
    rho_mt = (1 + c*rr + (5/3)*(rr.^2)) .* exp(-c*rr);
else
    coef   = 2^(1-nu) / gamma(nu);
    rho_mt = coef .* (rr.^nu) .* besselk(nu, rr);
end

% d>0: γ(d) = c0 + sigma2 * (1 - rho_nu)
g(mask) = c0 + s2 .* (1 - rho_mt);

end


% ---------- 其他小工具 ----------
function tf = in_lon_window_180(lon, lon0, halfwin_deg)
d = lon - lon0;
d = mod(d + 180, 360) - 180;  % wrap 到 [-180,180]
tf = abs(d) <= halfwin_deg;
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

% ----------------- 你现有的拟合函数（保持不变） -----------------
% 注意：fit 阶段可继续用 PSO_maternfit；但不要把函数句柄写进 MAT
function [ST2, Smodel, Tmodel] = fit_ST2_from_loc(ST_S, ST_T, ST_ST)
if istable(ST_T),  ST_T  = table2array(ST_T);  end
if istable(ST_S),  ST_S  = table2array(ST_S);  end
if istable(ST_ST), ST_ST = table2array(ST_ST); end

% —— 空间 —— %
if ~isempty(ST_S)
    s_lags_s = ST_S(:,1); exp_values_s = ST_S(:,2);
    [~, s_centers, s_grouped_values] = bin_grouping(s_lags_s, exp_values_s, 50);
    Smodel = PSO_maternfit(s_centers, s_grouped_values, 0.01, 0.8);
else
    Smodel = [];
end

% —— 时间（至少4个时间差） —— %
if ~isempty(ST_T)
    t_lags_t = ST_T(:,1); exp_values_t = ST_T(:,2);
    [t_unique, ~, t_group_idx] = unique(t_lags_t);
    t_grouped_values = accumarray(t_group_idx, exp_values_t, [], @mean);
    Tmodel = PSO_maternfit(t_unique, t_grouped_values, 0.01, 0.8);
else
    Tmodel = [];
end

% —— 时空 + α —— %
if isempty(ST_ST) || isempty(Smodel) || isempty(Tmodel)
    ST2 = []; return;
end

s_lags = ST_ST(:,1);
t_lags = ST_ST(:,2);
exp_values = ST_ST(:,3);

num_s_bins = 50;
[s_centers2, t_unique2, grouped_grid] = bin_grouping_2d(s_lags, t_lags, exp_values, num_s_bins);
[s_grid, t_grid] = meshgrid(s_centers2, t_unique2);

ys0 = Smodel.func(s_grid(:));    % 仅在拟合期用到函数句柄
y0t = Tmodel.func(t_grid(:));
Y   = grouped_grid(:);

valid = isfinite(ys0) & isfinite(y0t) & isfinite(Y);
ys0 = ys0(valid); y0t = y0t(valid); Y = Y(valid);

% Y = ys0 + y0t - α * (ys0 .* y0t) + ε
x = ys0 .* y0t;
r = Y - ys0 - y0t;
mask = (x > 0) & isfinite(r);

if ~any(mask)
    alpha_hat = 0;
else
    xm = x(mask); rm = r(mask);
    num = sum( xm .* (-rm) );
    den = sum( xm.^2 );
    alpha_hat = (den <= eps) * 0 + (den > eps) * (num/den);
end

sS = Smodel.sill; sT = Tmodel.sill;
alpha_hat = max(0, min(alpha_hat, min(1/max(sS,eps), 1/max(sT,eps))));

% 组装：注意这里不再设置 ST2.func（避免句柄被保存）
ST2.Smodel = Smodel;
ST2.Tmodel = Tmodel;
ST2.alpha  = alpha_hat;
ST2.nugget = Smodel.c0 + Tmodel.c0 - alpha_hat * Smodel.c0 * Tmodel.c0;
ST2.sill   = sS + sT - alpha_hat * sS * sT;
end
%% 大周期均值模型
function ST2 = fit_ST2_alpha_from_ST_only(ST_ST, Smodel, Tmodel)
% 仅根据 ST_ST 点 + 固定的 S/T 模型，闭式估计 alpha，并组装 ST2（不含 func）

if istable(ST_ST), ST_ST = table2array(ST_ST); end
if isempty(ST_ST)
    error('ST_ST is empty when fitting alpha.');
end

s_lags = ST_ST(:,1);
t_lags = ST_ST(:,2);
exp_values = ST_ST(:,3);

num_s_bins = 50;
[s_centers2, t_unique2, grouped_grid] = bin_grouping_2d(s_lags, t_lags, exp_values, num_s_bins);
[s_grid, t_grid] = meshgrid(s_centers2, t_unique2);

ys0 = Smodel.func(s_grid(:));     % 由块平均 S 模型给定
y0t = Tmodel.func(t_grid(:));     % 由块平均 T 模型给定
Y   = grouped_grid(:);

valid = isfinite(ys0) & isfinite(y0t) & isfinite(Y);
ys0 = ys0(valid); y0t = y0t(valid); Y = Y(valid);

x = ys0 .* y0t;
r = Y - ys0 - y0t;
mask = (x > 0) & isfinite(r);

if ~any(mask)
    alpha_hat = 0;
else
    xm = x(mask); rm = r(mask);
    num = sum( xm .* (-rm) );
    den = sum( xm.^2 );
    alpha_hat = (den <= eps) * 0 + (den > eps) * (num/den);
end

% 约束 α（与原版一致）
sS = Smodel.sill; sT = Tmodel.sill;
alpha_hat = max(0, min(alpha_hat, min(1/max(sS,eps), 1/max(sT,eps))));

% 组装 ST2（不设 func，避免句柄被存盘）
ST2.Smodel = Smodel;
ST2.Tmodel = Tmodel;
ST2.alpha  = alpha_hat;
ST2.nugget = Smodel.c0 + Tmodel.c0 - alpha_hat * Smodel.c0 * Tmodel.c0;
ST2.sill   = sS + sT - alpha_hat * sS * sT;
end
