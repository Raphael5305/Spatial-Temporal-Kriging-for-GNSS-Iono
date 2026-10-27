function S = variogram_spatial_fast(x, y, varargin)
% x = [lat, lon, epoch]; y = 值
% 输出 S = [distance_km, gamma, epoch]，按行堆叠所有 epoch 的成对结果

    % ---------------- 参数 ----------------
    prm.sth = 400;                 % 邻域半径 (km)
    prm.R   = 6378.137;            % 地球半径 (km)
    prm.use_single = false;        % 输出是否压成 single
    prm.max_pairs_per_epoch = inf; % 每个 epoch 限制最多对数（可用于下采样加速）
    prm = parseargs(prm, varargin{:});

    % -------------- 预处理 ---------------
    if size(y,1)~=size(x,1), error('x/y 行数不一致'); end
    ok = ~any(isnan(x),2) & ~isnan(y);
    x = x(ok,:); y = y(ok);

    epochs = unique(x(:,3));
    nE = numel(epochs);

    S_cells = cell(nE,1);

    % chord 半径（单位球上搜索半径），用来快速筛邻居
    theta_max = prm.sth / prm.R;                  % 弧度弧长
    r_chord   = 2*sin(0.5*theta_max);            % 单位球上 chord 距离

    for k = 1:nE
        ek = epochs(k);
        mask = x(:,3)==ek;
        if nnz(mask) < 2, continue; end

        lat = deg2rad(x(mask,1));
        lon = deg2rad(x(mask,2));
        yy  = y(mask);

        % ---- 构造单位球坐标（邻居搜索更快）----
        clat = cos(lat);
        X = [clat.*cos(lon), clat.*sin(lon), sin(lat)];  % [n×3] unit vectors

        % ---- 邻居搜索（优先用 rangesearch）----
        idxI = []; idxJ = [];
        if exist('rangesearch','file') == 2
            neigh = rangesearch(X, X, r_chord); % cell，每点的邻居（包含自己）
            % 只取上三角 i<j，避免重复
            n = numel(neigh);
            % 预估容量（保守）
            cap = 0;
            for i=1:n, cap = cap + max(0,numel(neigh{i})-1); end
            idxI = zeros(cap,1); idxJ = zeros(cap,1);
            c = 0;
            for i = 1:n
                js = neigh{i};
                js = js(js>i);
                m  = numel(js);
                if m>0
                    idxI(c+(1:m)) = i;
                    idxJ(c+(1:m)) = js(:);
                    c = c + m;
                end
            end
            if c<cap, idxI = idxI(1:c); idxJ = idxJ(1:c); end
        else
            % ---- 无 toolbox：向量化 fallback（块处理避免 OOM）----
            n = size(X,1);
            % 分块上三角
            blk = 5000;   % 视内存调整
            idxI = []; idxJ = [];
            for i1 = 1:blk:n
                i2 = min(n, i1+blk-1);
                % 仅与右上区域匹配
                j1 = i1; j2 = n;
                % 计算子块两两点的 chord 距离
                Dch = 1 - X(i1:i2,:)*X(j1:j2,:)';       % = (1 - dot), ∈[0,2]
                Dch = sqrt(max(0, 2*Dch));              % chord 距离
                % 小于阈值的对
                [ii,jj] = find(Dch <= r_chord + 1e-12);
                jj = jj + (j1-1);
                ii = ii + (i1-1);
                % 仅保留 i<j
                sel = ii < jj;
                idxI = [idxI; ii(sel)];
                idxJ = [idxJ; jj(sel)];
            end
        end

        if isempty(idxI), continue; end

        % ---- 精确大圆距离（向量化 haversine）----
        lat1 = lat(idxI); lon1 = lon(idxI);
        lat2 = lat(idxJ); lon2 = lon(idxJ);
        dkm  = great_circle_km(lat1,lon1,lat2,lon2, prm.R);

        % 只保留 d<=sth 的对（理论上已满足，这里再保险）
        keep = dkm <= prm.sth + 1e-9;
        if any(~keep)
            idxI = idxI(keep); idxJ = idxJ(keep);
            dkm  = dkm(keep);
        end
        if isempty(idxI), continue; end

        % ---- 计算 gamma = 0.5(y_i - y_j)^2 ----
        lam = 0.5*(yy(idxI) - yy(idxJ)).^2;

        % ---- 如需要，对过多的对数进行下采样 ----
        if numel(lam) > prm.max_pairs_per_epoch
            sel = randperm(numel(lam), prm.max_pairs_per_epoch);
            dkm = dkm(sel); lam = lam(sel);
        end

        if prm.use_single
            S_cells{k} = single([dkm, lam, repmat(single(ek), numel(lam),1)]);
        else
            S_cells{k} = [dkm, lam, repmat(ek, numel(lam),1)];
        end
    end

    if all(cellfun('isempty',S_cells))
        S = [];
    else
        S = vertcat(S_cells{:});
    end
end

% --------- 工具函数：向量化大圆距离（km）---------
function dkm = great_circle_km(lat1,lon1,lat2,lon2,R)
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    s = sin(0.5*dlat).^2 + cos(lat1).*cos(lat2).*sin(0.5*dlon).^2;
    ang = 2*atan2( sqrt(max(s,0)), sqrt(max(1-s,0)) );
    dkm = R * ang;
end
%% Subfunction: parseargs
function X = parseargs(X,varargin)

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if iscell(oldvalue)
            % Cell arrays must contain strings, and the new value must be
            % a string which appears in the list.
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings',fieldname));  %#ok
            end
            if ~ischar(newvalue)
                error(sprintf('New value for "%s" must be a string',fieldname));  %#ok
            end
            if isempty(find(strcmp(oldvalue,newvalue))) %#ok
                error(sprintf('"%s" is not allowed for field "%s"',newvalue,fieldname));  %#ok
            end
        elseif ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end

% Now find cell arrays which were not modified by the above process, and select
% the first string.
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X,fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue)
            error(sprintf('All allowed values for "%s" must be strings',fieldname)); %#ok
        elseif isempty(oldvalue)
            error(sprintf('Empty cell array not allowed for field "%s"',fieldname)); %#ok
        end
        X = setfield(X,fieldname,oldvalue{1}); %#ok
    end
end
end
