function exp_model = NLS_exfit(d, gamma_observed)
% NLS_exponentialfit 使用 fminsearchbnd 进行 Exponential 变差函数非线性最小二乘拟合
% 无权重、无绘图；返回字段与 NLS_maternfit 尽量保持一致以便平滑替换。
%
% 输入：
%   d              -- 滞后距离向量 (列向量/行向量均可)
%   gamma_observed -- 实验变差函数值
%
% 输出：
%   exp_model 结构体，字段包括：
%     .c0, .sigma2, .sill, .rho, .nu([]), .func
%     .obs, .gammahat, .residuals
%
% 依赖：
%   需确保 fminsearchbnd.m 位于路径中（与 Matern 版本一致）。

    % ---------- 输入预处理 ----------
    d = d(:);
    gamma_observed = gamma_observed(:);
    finite = isfinite(d) & isfinite(gamma_observed);
    d = d(finite);
    gamma_observed = gamma_observed(finite);

    if isempty(d) || isempty(gamma_observed) || var(gamma_observed) == 0
        exp_model = [];
        return;
    end

    % ---------- 降级优化分支（与 Matern 风格一致，给出闭式 sigma2 解） ----------
    % 规则：若 d 全为 30 的倍数 且 样本数 <= 4，则固定 c0=0, rho=180，仅估 sigma2（最小二乘闭式）
    tol = 1e-9;
    is_mult30 = all(abs(d/30 - round(d/30)) < tol);
    if is_mult30 && numel(d) <= 3
        c0_fix  = 0;
        rho_fix = 180;

        % phi(h) = 1 - exp(-h/rho)
        phi = zeros(size(d));
        idx = d > 0;
        phi(idx) = 1 - exp(-d(idx) / rho_fix);
        % 保持 PSO 风格：gamma(d==0) = 0；在 c0=0 情况下，phi(0)=0 符合实现

        nume = sum(phi .* (gamma_observed - c0_fix));
        deno = sum(phi.^2);
        if deno <= eps
            exp_model = [];
            return;
        end
        sigma2_hat = max(0, nume / deno);

        % 输出（字段名与上游保持一致）
        exp_model.c0     = c0_fix;        % nugget
        exp_model.sigma2 = sigma2_hat;    % partial sill
        exp_model.sill   = exp_model.c0 + exp_model.sigma2;
        exp_model.rho    = rho_fix;       % range/scale
        exp_model.nu     = [];            % 指数模型无 nu，置空以兼容
        exp_model.func   = @(dd) exponential_variogram_psostyle(dd, ...
                                 exp_model.c0, exp_model.sigma2, exp_model.rho);

        gammahat  = exp_model.func(d);
        residuals = gamma_observed - gammahat;
        exp_model.obs       = gamma_observed;
        exp_model.gammahat  = gammahat;
        exp_model.residuals = residuals;
        return;
    end

    % ---------- 常规 3 参数 NLS（fminsearchbnd） ----------
    if numel(d) < 4
        exp_model = [];
        return;
    end

    gmin = min(gamma_observed);
    gmax = max(gamma_observed);
    dmin = min(d);
    dmax = max(d);

    % 参数顺序: [c0, sigma2, rho]
    lb = [0,         0.1*(gmax-gmin), max(0.01*(dmax), 1e-6) ];
    ub = [gmax,      5*gmax,          max(5*dmax, 1500)        ];

    % 初始猜测（与 Matern 版本风格一致）
    sigma2_init = max(gmax - max(gmin,0), eps);
    rho_init    = 0.7 * max(d);
    c0_init     = max(gmin, 0);
    x0 = [c0_init, sigma2_init, rho_init];
    x0 = max(lb, min(ub, x0));

    % 目标函数（无权重 SSE）
    exp_fun = @(p, dd) exponential_variogram_psostyle(dd, p(1), p(2), p(3));
    obj_fun = @(p) sum( (exp_fun(p, d) - gamma_observed).^2 );

    options = optimset('Display','off', 'MaxFunEvals',3e6, 'MaxIter',3e6, ...
                       'TolFun',1e-6, 'TolX',1e-6);

    try
        [p_opt, ~, exitflag, ~] = fminsearchbnd(obj_fun, x0, lb, ub, options);
    catch
        warning('fminsearchbnd not found or failed, returning empty model.');
        exp_model = [];
        return;
    end

    if any(~isfinite(p_opt)) || exitflag <= 0
        exp_model = [];
        return;
    end

    % 输出（字段保持与你项目一致）
    exp_model.c0     = p_opt(1);                     % nugget
    exp_model.sigma2 = p_opt(2);                     % partial sill
    exp_model.sill   = p_opt(1) + p_opt(2);          % total sill
    exp_model.rho    = p_opt(3);                     % range/scale
    exp_model.nu     = [];                           % 指数模型无 nu
    exp_model.func   = @(dd) exponential_variogram_psostyle(dd, p_opt(1), p_opt(2), p_opt(3));

    gammahat  = exp_model.func(d);
    residuals = gamma_observed - gammahat;

    exp_model.obs       = gamma_observed;
    exp_model.gammahat  = gammahat;
    exp_model.residuals = residuals;
end


% ==== 指数变差函数（PSO 风格：d==0 → gamma(0)=0）====
function gamma_value = exponential_variogram_psostyle(d, c0, sigma2, rho)
    gamma_value = zeros(size(d));
    idx = d > 0;
    % 标准指数模型：gamma(h) = c0 + sigma2*(1 - exp(-h/rho))
    gamma_value(idx) = c0 + sigma2 .* (1 - exp(-d(idx) ./ rho));
    gamma_value(d == 0) = 0;  % 与你现有 PSO 实现保持一致
end
