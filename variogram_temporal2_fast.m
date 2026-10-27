function S = variogram_temporal2_fast(data, varargin)
% 快速版：在空间阈值内的点对上，计算不同历元差的时序半变差
% 输入  data = [lat, lon, epoch, y]
% 输出  S = [tau_seconds, gamma, prev_epoch, curr_epoch]
%
% 参数:
%   'sth'     (km)     空间阈值，默认 5
%   'n'       (steps)  时间窗口（按最小步长的 n 倍），默认 10
%   'R'       (km)     地球半径，默认 6378.137
%   'use_single'       是否把输出压成 single，默认 false
%   'max_pairs_per_ep' 每个 (prev,curr) 历元对最多抽样的点对数，默认 inf（不抽样）

% ----- 参数 -----
prm.sth = 5;
prm.n   = 10;
prm.R   = 6378.137;
prm.use_single = false;
prm.max_pairs_per_ep = inf;
prm = parseargs(prm, varargin{:});

% ----- 预处理 -----
ok = ~any(isnan(data),2);
data = data(ok,:);
if isempty(data), S = []; return; end

lat = data(:,1); lon = data(:,2); t = data(:,3); y = data(:,4);
epochs = unique(t);
if numel(epochs) < 2
    S = []; return;
end

% 推断时间步长（众数）
dt_unique = diff(sort(epochs));
if isempty(dt_unique)
    dt_step = 30; % 兜底
else
    dt_step = mode(round(dt_unique));
end
% 允许的窗口秒数
t_window = prm.n * dt_step;

% ----- 单位球坐标 & 半径邻居表 -----
% 邻居先用“弦长阈值”预筛，减少后续工作量
latr = deg2rad(lat); lonr = deg2rad(lon);
clat = cos(latr);
U = [clat.*cos(lonr), clat.*sin(lonr), sin(latr)];  % 单位向量

theta_max = prm.sth / prm.R;          % 弧长（弧度）
r_chord   = 2*sin(0.5*theta_max);     % 单位球 chord 距离阈值

if exist('rangesearch','file') == 2
    neigh = rangesearch(U, U, r_chord);  % cell: 每个点的邻居（含自己）
else
    % 无工具箱：用向量化块（内存可调）
    neigh = local_rangesearch_fallback(U, r_chord);
end

% ----- 建立 epoch → 索引 映射（在 data 的行号上） -----
[t_sorted, ord] = sort(t);
[t_unique, ~, ic] = unique(t_sorted);
groups = accumarray(ic, (1:numel(ic))', [], @(v){v});  % 每一唯一时刻在排序序列里的下标集合
% 为了快速从 epoch 值直接拿到原始行下标集合：
epoch_to_idx = containers.Map('KeyType','double','ValueType','any');
for i = 1:numel(t_unique)
    epoch_to_idx(t_unique(i)) = ord(groups{i});
end

% ----- 主循环：对每个当前历元，只匹配窗口内的过往历元 -----
Tbuf = {}; Lbuf = {}; Pbuf = {}; Cbuf = {};
for k = 2:numel(epochs)
    curr = epochs(k);
    A = epoch_to_idx(curr);
    if isempty(A), continue; end

    % 窗口内的过往历元（>= curr - t_window 且 < curr）
    lo = curr - t_window;
    % epochs 已排序，二分找起点
    m_start = find(epochs >= lo, 1, 'first');
    if isempty(m_start), m_start = 1; end

    for m = m_start:k-1
        prev = epochs(m);
        B = epoch_to_idx(prev);
        if isempty(B), continue; end

        % ---- 把 B 做成逻辑表，快速过滤近邻 ----
        isB = false(numel(y),1); isB(B) = true;

        % A: 当前历元索引集合；B: 过往历元索引集合
        % neigh: 全局邻接表（cell，含自身）；isB: 属于 B 的逻辑向量

        % 1) 先得到 A 中每个点在 B 中的邻居（列向量），并剔除空单元
        pair_j = cellfun(@(lst) reshape(lst(isB(lst)), [], 1), neigh(A), 'uni', false);
        nz = ~cellfun('isempty', pair_j);
        pair_j = pair_j(nz);
        A_nz   = A(nz);

        % 2) 为每个单元生成对应数量的 i 索引（列向量）
        pair_i = cellfun(@(i,js) repmat(i, numel(js), 1), num2cell(A_nz), pair_j, 'uni', false);

        % 3) 拼接成完整索引对
        I = vertcat(pair_i{:});
        J = vertcat(pair_j{:});
        if isempty(I)
            % 这个 prev-curr 历元对没有任何空间近邻
            continue;
        end

        % 4) —— 严格用大圆距离判定 ≤ sth（与原暴力版一致）——
        % 向量化 haversine（与 distance('gc')*R 数值等价）
        dlat = latr(J) - latr(I);
        dlon = lonr(J) - lonr(I);
        s = sin(0.5*dlat).^2 + cos(latr(I)).*cos(latr(J)).*sin(0.5*dlon).^2;
        ang = 2*atan2( sqrt(max(s,0)), sqrt(max(1-s,0)) );  % 弧度
        dkm = prm.R * ang;                                      % 公里

        keep = dkm <= prm.sth + 1e-12;   % 小容差确保边界一致
        if ~any(keep)
            continue;
        end
        I = I(keep);
        J = J(keep);

        % ---- 半变差 ----
        lam = 0.5*(y(I) - y(J)).^2;

        % ---- 记录 ----
        tau = curr - prev;
        Tbuf{end+1} = repmat(tau, numel(lam), 1);
        Lbuf{end+1} = lam;
        Pbuf{end+1} = repmat(prev, numel(lam), 1);
        Cbuf{end+1} = repmat(curr, numel(lam), 1);
    end
end

if isempty(Lbuf)
    S = [];
else
    TemporalDistance = vertcat(Tbuf{:});
    Variogram        = vertcat(Lbuf{:});
    PreviousEpoch    = vertcat(Pbuf{:});
    CurrentEpoch     = vertcat(Cbuf{:});
    S = [TemporalDistance, Variogram, PreviousEpoch, CurrentEpoch];
    if prm.use_single, S = single(S); end
end
end

% ===== 工具函数：无 Statistics Toolbox 时的 rangesearch 替代 =====
function neigh = local_rangesearch_fallback(U, r_chord)
n = size(U,1);
neigh = cell(n,1);
% 块大小可按内存调大（越大越快，内存占用增加）
blk = 8000;
for i1 = 1:blk:n
    i2 = min(n, i1+blk-1);
    % 与整个集合做点积（两层块可进一步降低峰值内存，这里先简单实现）
    D = 1 - U(i1:i2,:)*U.';         % 1 - dot ∈ [0,2]
    D = sqrt(max(0, 2*D));          % chord 距离
    for ii = 1:(i2-i1+1)
        js = find(D(ii,:) <= r_chord + 1e-12);
        neigh{i1+ii-1} = js(:).';
    end
end
end

% ===== 解析参数（保持你原来的 parseargs 语义） =====
function X = parseargs(X, varargin)
remaining = nargin - 1;
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
while remaining >= 2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname, fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X, fieldname); %#ok
        newvalue = varargin{count + 1};
        if iscell(oldvalue)
            if ~iscellstr(oldvalue), error('All allowed values for "%s" must be strings', fieldname); end
            if ~ischar(newvalue),   error('New value for "%s" must be a string', fieldname); end
            if isempty(find(strcmp(oldvalue, newvalue),1))
                error('"%s" is not allowed for field "%s"', newvalue, fieldname);
            end
        elseif ~isempty(oldvalue)
            if ~strcmp(class(oldvalue), class(newvalue))
                error('Cannot change class of field "%s" from "%s" to "%s"', fieldname, class(oldvalue), class(newvalue));
            elseif numel(oldvalue)==1 && numel(newvalue)~=1
                error('New value for "%s" must be a scalar', fieldname);
            end
        end
        X = setfield(X, fieldname, newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
if remaining ~= 0, error('Odd number of arguments supplied. Name-value pairs required'); end
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X, fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue), error('All allowed values for "%s" must be strings', fieldname); end
        if isempty(oldvalue),    error('Empty cell array not allowed for field "%s"', fieldname); end
        X = setfield(X, fieldname, oldvalue{1}); %#ok
    end
end
end
