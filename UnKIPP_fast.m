function results = UnKIPP_fast(R, N, IPPVD, samp_int, dq)
% simplethreat_fast — 无 data deprivation 版
% 行步长选 target（与原版一致），同一历元内按距离≤R_max 取样；
% 自适应确定 Rfitt（确保样本数≥N_min），用默认 UnK 模型做一次估计。
%
% 输出(每行 18 列):
% [ epoch, sat, lat, lon, y_true, zi, chi, chi_irreg, R_irreg, sigma_IPP, sigma_inflated, ...
%   Rfitt, RCM1, R_T, TCM, NaN, NaN, NaN ]

if nargin < 5, dq = []; end %#ok<NASGU>

% ==== default UnK model（与原版一致）====
defaultUnK.func        = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
defaultUnK.sill        = 0.85;
defaultUnK.model       = 'exponential';
defaultUnK.nugget      = 0.05;
defaultUnK.partialsill = 0.8;
defaultUnK.range       = 8000;   % km
defaultfit             = [defaultUnK.range, defaultUnK.partialsill, defaultUnK.nugget];

R_min = R(1); R_max = R(2);
N_min = N(1); N_target = N(2); %#ok<NASGU>

wgs84   = wgs84Ellipsoid('kilometer');
height  = 350;          % km
R_earth = 6378.137;     % km
DEG2RAD = pi/180;

PF    = IPPVD;          % [epoch, sat, lat, lon, ..., y]
nrows = size(PF,1);

% === 仍构建 epoch→行号索引，但目标点采样按“行步长” ===
[~, ~, grp] = unique(PF(:,1), 'stable');        % grp(i) = epoch id of row i
idx_cells   = accumarray(grp, (1:nrows).', [], @(v){v});  % 每个epoch的所有行号

% —— 采样对象：严格等同原版 ——
target_rows = 1:samp_int:nrows;
data_est    = PF(target_rows, :);
if isempty(data_est)
    results = zeros(0,18);
    return;
end

results = zeros(numel(target_rows), 18);
save_row = 0;

for t = 1:numel(target_rows)
    irow     = target_rows(t);
    est_time = PF(irow, 1);

    % 同一历元的所有观测
    gid        = grp(irow);
    rows_epoch = idx_cells{gid};
    if isempty(rows_epoch), continue; end

    % 目标坐标
    lat1 = PF(irow,3);
    lon1 = PF(irow,4);
    if lat1<=22.5||lat1>=52.5||lon1>=-62.5||lon1<=-137.5
        continue;
    end

    % 候选坐标
    lat2 = PF(rows_epoch,3);
    lon2 = PF(rows_epoch,4);

    % 大圆距离
    dlat = (lat2 - lat1)*DEG2RAD;
    dlon = (lon2 - lon1)*DEG2RAD;
    a = sin(0.5*dlat).^2 + cos(lat1*DEG2RAD) .* cos(lat2*DEG2RAD) .* sin(0.5*dlon).^2;
    Dist_full = R_earth .* (2*atan2(sqrt(a), sqrt(1 - a)));

    % 按距离升序（去掉自身）
    [Dist_sorted, ord] = sort(Dist_full, 'ascend');
    rows_sorted = rows_epoch(ord);
    if ~isempty(Dist_sorted)
        Dist_sorted(1) = [];      % 删除与 target 重合点
        rows_sorted(1) = [];
    end

  % WAAS 四步规则确定 Rfitt
        nRmin = sum(Dist_sorted <= R_min);
        if nRmin >= N_target
            Rfitt = R_min;
        elseif numel(Dist_sorted) >= N_target && Dist_sorted(N_target) <= R_max
            Rfitt = Dist_sorted(N_target);
        else
            Rfitt = R_max;
            if sum(Dist_sorted <= Rfitt) < N_min
                continue; % Rmax 内样本不足
            end
        end

    sel = (Dist_sorted <= Rfitt);
    if ~any(sel), continue; end
    id_sel = rows_sorted(sel);
    PF_samples = PF(id_sel, :);
    if size(PF_samples,1) < N_min, continue; end

    % East/North（相对 target 位置）
    East2  = zeros(size(PF_samples,1),1);
    North2 = zeros(size(PF_samples,1),1);
    for k2 = 1:size(PF_samples,1)
        [E2, N2, ~] = geodetic2enu(PF_samples(k2,3), PF_samples(k2,4), height, ...
                                   lat1, lon1, height, wgs84, 'degrees');
        East2(k2)  = E2;  North2(k2) = N2;
    end

    % 距离矩阵（UnK：无时间维）
    ns = size(PF_samples,1);
    dist_n0 = distance('gc', PF_samples(:,3), PF_samples(:,4), lat1, lon1, R_earth);
    [Mg, Ng] = meshgrid(1:ns);
    dist_nn = reshape(distance('gc', ...
        PF_samples(Mg(:),3), PF_samples(Mg(:),4), PF_samples(Ng(:),3), PF_samples(Ng(:),4), R_earth), ns, ns);

    % RCM/时间量（单历元：R_T=0, TCM=0）
    Rcentroid1 = mean(dist_n0);
    RCM1 = Rcentroid1 / max(Rfitt, eps);
    R_T  = 0; TCM = 0;

    % UnK 估计（与 NLSSTGIVx2 的 UnK_threat 输出对应）
    data_values = PF_samples(:,5);
    [X0, Y0] = meshgrid(0);
    UnK_result = Unk_threat(defaultUnK, East2, North2, data_values, X0, Y0, defaultfit, dist_nn, dist_n0);
    % UnK_result: [zi, chi, chi_irreg, R_irreg, sigma_IPP, sigma_inflated]

    % 写结果：保持 18 列；末尾 pad3 预留占位方便兼容
    target_point = PF(irow, 1:5);  % [epoch, sat, lat, lon, y_true]
    pad3 = [NaN, NaN, NaN];

    save_row = save_row + 1;
    if save_row > size(results,1)
        results = [results; zeros(1024, 18)]; %#ok<AGROW>
    end
    results(save_row, :) = [ target_point, UnK_result(:).', Rfitt, RCM1, R_T, TCM, pad3 ];

    % 进度
    if ~isempty(dq) && mod(t, 10000) == 0
        task = getCurrentTask(); if ~isempty(task), tid = task.ID; else, tid = 0; end
        send(dq, sprintf('线程 %d | target %d/%d (row %d)', tid, t, numel(target_rows), irow));
    end
end

% 裁剪
results = results(1:save_row, :);
end
