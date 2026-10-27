function matern_model=PSO_maternfit(d,gamma_observed,mse_th,Rs_th)
% Explanation: d-variogram factors
% observed-the calculated variogram values
% mse_th-the threshold of mse error
% Rs_th-the threshold of Rsquare
%数据不足（尤其是time），不收敛，解算失败时返回空数组

    % 列向量化 + 有限值过滤
    d = d(:);
    gamma_observed = gamma_observed(:);
    finite = isfinite(d) & isfinite(gamma_observed);
    d = d(finite);
    gamma_observed = gamma_observed(finite);

    % 样本与滞后值丰富度检查（不可辨识时直接返回空）
    if numel(d) < 5 || numel(unique(round(d, 6))) < 5
        matern_model = [];
        return;
    end
   % 观测方差为 0，无法拟合
    if var(gamma_observed) == 0
        matern_model = [];
        return;
    end
    % 使用 PSO 拟合 Matern 变差函数
    lb = [0, 0.1*(max(gamma_observed)-min(gamma_observed)), 0.1*(max(d)-min(d)), 0.01];  % 参数下界 [c0, sigma2, rho, nu]
    ub = [max(gamma_observed), 5*max(gamma_observed), max(5*max(d),1500), 20.0];  % 参数上界 [c0, sigma2, rho, nu]

    
    options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 4000,'MaxStallIterations', 100, 'FunctionTolerance', 1e-8,'Display','off');
    %% 如果考虑非线性最小二乘初始化PSO参数加快收敛，可以使用以下，需要修改： % 初始参数启发式猜测
    % sigma2_init = max(gamma_observed);
    % idx = find(gamma_observed >= 0.7*sigma2_init, 1);
    % if isempty(idx)
    %     rho_init = max(d)/2;
    % else
    %     rho_init = d(idx);
    % end
    % nu_init = 1.5; % 常用初始值
    % 
    % % 非线性最小二乘拟合
    % lsq_lb = [0.1, 0.1, 0.1];
    % lsq_ub = [10, 10, 10];
    % params_init = [sigma2_init, rho_init, nu_init];
    % matern_fun = @(params, d) matern_variogram(d, params(1), params(2), params(3));
    % 
    % try
    %     options_lsq = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 100);
    %     params_lsq = lsqcurvefit(matern_fun, params_init, d, gamma_observed, lsq_lb, lsq_ub, options_lsq);
    % catch
    %     params_lsq = params_init; % 失败时使用初始猜测
    % end
    % 
    % % 调整PSO的上下界至非线性结果附近
    % range_scale = 0.5; % 搜索范围扩展50%
    % lb = max([0.1, 0.1, 0.1], params_lsq * (1 - range_scale));
    % ub = min([10.0, 10.0, 10.0], params_lsq * (1 + range_scale));
    % 
    % % 设置PSO初始种群（包含非线性结果）
    % swarmSize = 30;
    % initialSwarm = repmat(params_lsq, ceil(swarmSize/3), 1); % 1/3粒子来自非线性结果
    % initialSwarm = [initialSwarm; rand(swarmSize - size(initialSwarm,1), 3) .* (ub - lb) + lb];
    % 
    % options = optimoptions('particleswarm', ...
    %                       'SwarmSize', swarmSize, ...
    %                       'MaxIterations', 200, ...
    %                       'InitialSwarmMatrix', initialSwarm, ...
    %                       'Display', 'off');
    %%  
    try
        [params_opt, mse_opt,exitflag] = particleswarm(@(params) objective_function(params, d, gamma_observed), 4, lb, ub, options);
    catch
        matern_model=[];
        return;
    end

    % 参数合法性与收敛性判定
    if any(~isfinite(params_opt)) || ~isfinite(mse_opt) || exitflag <= 0
        matern_model = [];
        return;
    end
    % 将拟合好的参数封装为模型
    matern_model.c0 = params_opt(1);      % nugget effect
    matern_model.sigma2 = params_opt(2);  % partialsill
    matern_model.sill=params_opt(1)+params_opt(2);
    matern_model.rho = params_opt(3);
    matern_model.nu = params_opt(4);
    matern_model.func = @(dd) matern_variogram(dd, matern_model.c0, matern_model.sigma2, matern_model.rho, matern_model.nu);
    % 计算预测值、残差、协方差和 R 平方
    gammahat = matern_model.func(d); % 预测值
    residuals = gamma_observed - gammahat; % 残差

    % R^2 计算健壮化
    vhat = var(gammahat);
    vobs = var(gamma_observed);
    if vhat == 0 || vobs == 0
        Rsq = NaN;
        C   = [NaN NaN; NaN NaN];
    else
        C   = cov(gammahat, gamma_observed);
        Rsq = (C(1,2)^2) / (vhat * vobs);
    end
    
    COVyhaty = C; % 协方差
    Rsquare = Rsq;% R 平方
    % 保存详细结果到结构体
    matern_model.obs = gamma_observed; % 实验值
    matern_model.gammahat = gammahat; % 预测值
    matern_model.residuals = residuals; % 残差
    matern_model.COVyhaty = COVyhaty; % 协方差
    matern_model.Rsquare = Rsquare; % R 平方

    % if mse_opt>=mse_th || Rsquare<=Rs_th
    %     disp(['mse:' num2str(mse_opt) 'threshold:' num2str(mse_th)]);
    %     disp(['Rs:' num2str(Rsquare) 'threshold:' num2str(Rs_th)]);
    % end
end

% Matern 变差函数
function gamma_value = matern_variogram(d,c0,sigma2, rho, nu)
    gamma_value = zeros(size(d));
    idx = d > 0;
    term1 = (d(idx) / rho).^nu;
    term2 = besselk(nu, d(idx) / rho);
    gamma_value(idx) = c0 + sigma2 * (1 - (2^(1 - nu) / gamma(nu)) * term1 .* term2);
    gamma_value(d == 0) = 0;  % 确保 d=0 时 gamma=0
end

% 目标函数（均方误差）
function mse = objective_function(params, d, gamma_observed)
    c0 = params(1);      % nugget effect
    sigma2 = params(2);  % partial sill
    rho = params(3);     % range
    nu = params(4);      % smoothness
    gamma_predicted = matern_variogram(d, c0, sigma2, rho, nu);
    mse = mean((gamma_predicted - gamma_observed).^2);
end