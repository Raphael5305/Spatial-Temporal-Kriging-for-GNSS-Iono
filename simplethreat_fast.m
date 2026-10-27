function results = simplethreat_fast(R, N, IPPVD, samp_int, Rin_int, modelType, mode, selfmodel, dq)
% simplethreat_fast
% 并行版（parfor），保持与原版结果和顺序一致
%
% 新增输入：
% modelType:
%   'gaussian' 或 'exponential'
%
% mode:
%   'default'  使用内置默认模型参数
%   'self'     使用 selfmodel
%
% selfmodel 结构体要求：
%   selfmodel.c0      -> nugget
%   selfmodel.sigma2  -> partial sill
%   selfmodel.sill    -> sill
%   selfmodel.rho     -> range

if nargin < 6 || isempty(modelType)
    modelType = 'exponential';
end
if nargin < 7 || isempty(mode)
    mode = 'default';
end
if nargin < 8
    selfmodel = [];
end
if nargin < 9
    dq = [];
end

modelType = lower(string(modelType));
mode      = lower(string(mode));

%% ==== build UnK model ====
switch mode
    case "default"
        % 默认参数
        UnKmodel.sill        = 1.25 * 1.25;
        UnKmodel.nugget      = 0.45 * 0.45;
        UnKmodel.partialsill = UnKmodel.sill - UnKmodel.nugget;
        UnKmodel.range       = 8000;   % km

    case "self"
        if isempty(selfmodel)
            error('When mode = ''self'', input selfmodel is required.');
        end

        requiredFields = {'c0','sigma2','sill','rho'};
        for k = 1:numel(requiredFields)
            if ~isfield(selfmodel, requiredFields{k})
                error('selfmodel must contain fields: c0, sigma2, sill, rho.');
            end
        end

        UnKmodel.nugget      = selfmodel.c0;
        UnKmodel.partialsill = selfmodel.sigma2;
        UnKmodel.sill        = selfmodel.sill;
        UnKmodel.range       = selfmodel.rho;

        tol = 1e-10;
        if abs(UnKmodel.sill - (UnKmodel.nugget + UnKmodel.partialsill)) > tol
            warning('selfmodel.sill ~= selfmodel.c0 + selfmodel.sigma2. The actual kriging uses nugget and partial sill.');
        end

    otherwise
        error('mode must be ''default'' or ''self''.');
end

%% ==== variogram type ====
switch modelType
    case "exponential"
        UnKmodel.func  = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
        UnKmodel.model = 'exponential';

    case "gaussian"
        UnKmodel.func  = @(b,h) b(2) * (1 - exp(-(h ./ b(1)).^2));
        UnKmodel.model = 'gaussian';

    otherwise
        error('modelType must be ''gaussian'' or ''exponential''.');
end

fitpara = [UnKmodel.range, UnKmodel.partialsill, UnKmodel.nugget];

R_min = R(1); 
R_max = R(2);
N_min = N(1); %#ok<NASGU>
N_target = N(2); %#ok<NASGU>

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350;          % km
R_earth = 6378.137;     % km
DEG2RAD = pi / 180;

PF    = IPPVD;          % [epoch, sat, lat, lon, ..., y]
nrows = size(PF,1);

% === 构建 epoch→行号索引 ===
[~, ~, grp] = unique(PF(:,1), 'stable');
idx_cells   = accumarray(grp, (1:nrows).', [], @(v){v});

% === 采样对象 ===
target_rows = 1:samp_int:nrows;
data_est    = PF(target_rows, :); %#ok<NASGU>
if isempty(target_rows)
    results = zeros(0,19);
    return;
end

% === 结果容量 ===
R_in_values = 0:Rin_int:2100;
Rw_values   = (0.1:0.2:0.9) * 2100;
per_target_cap = numel(R_in_values) * numel(Rw_values);

results_cell = cell(numel(target_rows), 1);

%% ===== 并行主循环 =====
for t = 1:numel(target_rows)
    irow     = target_rows(t);
    est_time = PF(irow, 1); %#ok<NASGU>

    gid        = grp(irow);
    rows_epoch = idx_cells{gid};
    if isempty(rows_epoch)
        results_cell{t} = zeros(0,19);
        continue;
    end

    lat1 = PF(irow,3);
    lon1 = PF(irow,4);
    cos_lat1 = cos(lat1 * DEG2RAD);

    lat2 = PF(rows_epoch,3);
    lon2 = PF(rows_epoch,4);

    % 向量化大圆距离
    dlat = (lat2 - lat1) * DEG2RAD;
    dlon = (lon2 - lon1) * DEG2RAD;
    a = sin(0.5*dlat).^2 + cos_lat1 .* cos(lat2*DEG2RAD) .* sin(0.5*dlon).^2;
    Dist_full = R_earth .* (2 * atan2(sqrt(a), sqrt(1 - a)));

    % 排序并去掉自点
    [Dist_sorted, ord] = sort(Dist_full, 'ascend');
    rows_sorted = rows_epoch(ord);

    if ~isempty(Dist_sorted)
        Dist_sorted(1) = [];
        rows_sorted(1) = [];
    end

    sel = (Dist_sorted <= R_max);
    if ~any(sel)
        results_cell{t} = zeros(0,20);
        continue;
    end

    idsearch_initial = rows_sorted(sel);
    PF_samples = PF(idsearch_initial, :);
    target_point = PF(irow, 1:5);

    buf  = zeros(per_target_cap, 19);
    used = 0;

    for r = 1:numel(R_in_values)
        R_in = R_in_values(r);

        for w = 1:numel(Rw_values)
            R_out = R_in + Rw_values(w);

            result_once = depriveAndEstimate2( ...
                R, N, R_in, R_out, ...
                PF_samples, target_point, ...
                UnKmodel, fitpara, wgs84, height, R_earth);

            if ~all(isnan(result_once))
                used = used + 1;
                buf(used, :) = [target_point, result_once];
            end
        end
    end

    if used > 0
        results_cell{t} = buf(1:used, :);
    else
        results_cell{t} = zeros(0,19);
    end

    if ~isempty(dq) && mod(t, 300) == 0
        task = getCurrentTask();
        if ~isempty(task)
            tid = task.ID;
        else
            tid = 0;
        end
        send(dq, sprintf('UnK线程 %d | target %d/%d (row %d)', tid, t, numel(target_rows), irow));
    end
end

results = vertcat(results_cell{:});

end