function matern_model = NLS_maternfit(d, gamma_observed)
% NLS_maternfit 使用 fminsearchbnd 进行 Matern 变差函数非线性最小二乘拟合
% 无权重、无绘图，输出结构与当前版本一致。
%
% 输入：
%   d              -- 滞后距离向量
%   gamma_observed -- 实验变差函数值
%
% 输出：
%   matern_model 结构体，字段包括：
%     .c0, .sigma2, .sill, .rho, .nu, .func
%     .obs, .gammahat, .residuals
%
% 依赖：
%   需确保 fminsearchbnd.m 位于路径中。

    % ---------- 输入预处理 ----------
    d = d(:);
    gamma_observed = gamma_observed(:);
    finite = isfinite(d) & isfinite(gamma_observed);
    d = d(finite);
    gamma_observed = gamma_observed(finite);

    if isempty(d) || isempty(gamma_observed) || var(gamma_observed)==0
        matern_model = [];
        return;
    end

    % ---------- 降级优化分支 ----------
    % 规则：若 d 全为 30 的倍数 且 样本数 <= 3，则固定 nu=0.5, c0=0, rho=60，仅估 sigma2（闭式解）
    tol = 1e-9;
    is_mult30 = all( abs(d/30 - round(d/30)) < tol );
    if is_mult30 && numel(d) <= 4
        c0_fix  = 0;
        nu_fix  = 0.5;
        rho_fix = 180;

        % phi(h) = 1 - (2^(1-nu)/Gamma(nu)) * (h/rho)^nu * K_nu(h/rho)
        phi = zeros(size(d));
        idx = d > 0;
        if any(idx)
            term1 = (d(idx)/rho_fix).^nu_fix;
            term2 = besselk(nu_fix, d(idx)/rho_fix);
            phi(idx) = 1 - (2^(1 - nu_fix)/gamma(nu_fix)) * term1 .* term2;
        end
        % 与 PSO 一致：gamma(d==0) = 0；在 c0=0 情况下，phi(0) 取 0 符合该实现

        nume = sum( phi .* (gamma_observed - c0_fix) );
        deno = sum( phi.^2 );
        if deno <= eps
            matern_model = [];
            return;
        end
        sigma2_hat = max(0, nume / deno);

        % 输出（与当前版本字段一致）
        matern_model.c0     = c0_fix;                 % nugget
        matern_model.sigma2 = sigma2_hat;             % partial sill
        matern_model.sill   = matern_model.c0 + matern_model.sigma2;
        matern_model.rho    = rho_fix;                % range/scale
        matern_model.nu     = nu_fix;                 % smoothness
        matern_model.func   = @(dd) matern_variogram_psostyle(dd, ...
                                   matern_model.c0, matern_model.sigma2, ...
                                   matern_model.rho, matern_model.nu);

        gammahat  = matern_model.func(d);
        residuals = gamma_observed - gammahat;
        matern_model.obs       = gamma_observed;
        matern_model.gammahat  = gammahat;
        matern_model.residuals = residuals;
        return;
    end

    % ---------- 常规 4 参数 NLS（fminsearchbnd） ----------
    if numel(d) < 5
        matern_model = [];
        return;
    end

    gmin = min(gamma_observed);
    gmax = max(gamma_observed);
    dmin = min(d);
    dmax = max(d);

    lb = [0, 0.1*(gmax-gmin), max(0.1*(dmax-dmin)), 0.1];  % [c0, sigma2, rho, nu]
    ub = [gmax, 5*gmax, max(5*dmax, 1500), 20.0];

    % 初始猜测（保持你当前版本设定）
    sigma2_init = gmax;
    rho_init    = 0.7*max(d);
    nu_init     = 0.5;
    c0_init     = max(gmin, 0);
    x0 = [c0_init, sigma2_init, rho_init, nu_init];
    x0 = max(lb, min(ub, x0));

    % 目标函数（无权重 SSE）
    matern_fun = @(p, dd) matern_variogram_psostyle(dd, p(1), p(2), p(3), p(4));
    obj_fun    = @(p) sum((matern_fun(p, d) - gamma_observed).^2);

    options = optimset('Display','off', 'MaxFunEvals',3e6, 'MaxIter',3e6, ...
                       'TolFun',1e-6,'TolX',1e-6);

    try
        [p_opt, ~, exitflag, ~] = fminsearchbnd(obj_fun, x0, lb, ub, options);
    catch
        warning('fminsearchbnd not found or failed, returning empty model.');
        matern_model = [];
        return;
    end

    if any(~isfinite(p_opt)) || exitflag <= 0
        matern_model = [];
        return;
    end

    % 输出（与当前版本字段一致）
    matern_model.c0     = p_opt(1);                 % nugget
    matern_model.sigma2 = p_opt(2);                 % partial sill
    matern_model.sill   = p_opt(1) + p_opt(2);      % total sill
    matern_model.rho    = p_opt(3);                 % range/scale
    matern_model.nu     = p_opt(4);                 % smoothness
    matern_model.func   = @(dd) matern_variogram_psostyle(dd, p_opt(1), p_opt(2), p_opt(3), p_opt(4));

    gammahat  = matern_model.func(d);
    residuals = gamma_observed - gammahat;

    matern_model.obs       = gamma_observed;
    matern_model.gammahat  = gammahat;
    matern_model.residuals = residuals;
end


% ==== 与 PSO 版本一致的 Matern 变差函数（保持 d==0 → 0） ====
function gamma_value = matern_variogram_psostyle(d, c0, sigma2, rho, nu)
    gamma_value = zeros(size(d));
    idx = d > 0;
    term1 = (d(idx) / rho).^nu;
    term2 = besselk(nu, d(idx) / rho);
    gamma_value(idx) = c0 + sigma2 * (1 - (2^(1 - nu) / gamma(nu)) * term1 .* term2);
    gamma_value(d == 0) = 0;  % 与 PSO_maternfit 保持一致
end
