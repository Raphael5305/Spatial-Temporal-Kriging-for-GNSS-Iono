function S = variogram_ST2_fast(x, y, varargin)
% 快速/等价版时空变差：
% 输入 x=[lat,lon,epoch]，y=值
% 输出 S=[distance_km, tau_s, gamma, prev_epoch, curr_epoch]

    % ---------- 参数 ----------
    prm.sth = 400;
    prm.n   = 5;                   % 30*prm.n 秒窗口
    prm.R   = 6378.137;
    prm.use_single = false;
    prm.max_pairs_per_epair = inf;
    prm = parseargs(prm, varargin{:});

    % ---------- 预处理 ----------
    if size(y,1) ~= size(x,1), error('x/y 行数不一致'); end
    ok = ~any(isnan(x),2) & ~isnan(y);
    x = x(ok,:); y = y(ok);
    if isempty(x), S=[]; return; end

    epochs = unique(x(:,3));
    epochs = sort(epochs);
    nE = numel(epochs);
    max_dt = 30 * prm.n;

    % ---------- 单位球坐标 & 空间近邻（半径=sth km） ----------
    lat  = x(:,1);  lon = x(:,2);  t = x(:,3);
    latr = deg2rad(lat); lonr = deg2rad(lon);
    clat = cos(latr);
    U = [clat.*cos(lonr), clat.*sin(lonr), sin(latr)];  % 单位球向量

    theta_max = prm.sth / prm.R;          % 弧度
    r_chord   = 2*sin(0.5*theta_max);     % 单位球弦长阈值（仅预筛）

    if exist('rangesearch','file') == 2
        neigh = rangesearch(U, U, r_chord);  % cell: 每点近邻（含自身）
    else
        neigh = local_rangesearch_fallback(U, r_chord);
    end

    % ---------- 历元 → 行索引 映射 ----------
    idx_by_epoch = cell(nE,1);
    for k=1:nE
        idx_by_epoch{k} = find(t==epochs(k));
    end

    % ---------- 结果缓冲 ----------
    Dbuf = {};  Tbuf = {};  Lbuf = {};  Pbuf = {};  Cbuf = {};

    % ========== (A) 同历元 Δt=0 ==========
    for k = 1:nE
        Ik = idx_by_epoch{k};
        if numel(Ik) < 2, continue; end

        isIk = false(numel(y),1); isIk(Ik) = true;

        % A 内近邻（列向量，剔空）
        pair_j = cellfun(@(lst) reshape(lst(isIk(lst)), [], 1), neigh(Ik), 'uni', false);
        nz = ~cellfun('isempty', pair_j);
        if ~any(nz), continue; end
        pair_j = pair_j(nz);
        Ai     = Ik(nz);

        % 展开为成对索引（列向量）
        pair_i = cellfun(@(i,js) repmat(i, numel(js), 1), num2cell(Ai), pair_j, 'uni', false);
        I = vcat_cols(pair_i);    % 安全拼接
        J = vcat_cols(pair_j);
        if isempty(I), continue; end

        % 上三角 i<j
        keep = I < J;
        I = I(keep); J = J(keep);
        if isempty(I), continue; end

        % 严格大圆距离筛 ≤ sth
        dkm = great_circle_km(latr(I), lonr(I), latr(J), lonr(J), prm.R);
        keep = dkm <= prm.sth + 1e-12;
        if ~any(keep), continue; end
        I = I(keep); J = J(keep); dkm = dkm(keep);

        lam = 0.5*(y(I) - y(J)).^2;
        ek = epochs(k);  % ← 修正为 MATLAB 索引
        Dbuf{end+1} = as_col(dkm);
        Tbuf{end+1} = zeros(numel(lam),1,'like',lam);
        Lbuf{end+1} = as_col(lam);
        Pbuf{end+1} = repmat(ek, numel(lam),1);
        Cbuf{end+1} = repmat(ek, numel(lam),1);
    end

    % ========== (B) 跨历元 0<Δt≤max_dt ==========
    for k = 2:nE
        t_cur = epochs(k);
        Ik = idx_by_epoch{k};
        if isempty(Ik), continue; end

        m_start = find(epochs >= t_cur - max_dt, 1, 'first');
        if isempty(m_start), m_start = 1; end

        for m = m_start:(k-1)
            t_prev = epochs(m);
            Im = idx_by_epoch{m};
            if isempty(Im), continue; end

            isB = false(numel(y),1); isB(Im) = true;

            pair_j = cellfun(@(lst) reshape(lst(isB(lst)), [], 1), neigh(Ik), 'uni', false);
            nz = ~cellfun('isempty', pair_j);
            if ~any(nz), continue; end
            pair_j = pair_j(nz);
            Ai     = Ik(nz);

            pair_i = cellfun(@(i,js) repmat(i, numel(js), 1), num2cell(Ai), pair_j, 'uni', false);
            I = vcat_cols(pair_i);
            J = vcat_cols(pair_j);
            if isempty(I), continue; end

            % 抽样（如有需要，且注意 0 对的情况）
            if numel(I) > prm.max_pairs_per_epair
                sel = randperm(numel(I), prm.max_pairs_per_epair);
                I = I(sel); J = J(sel);
            end
            if isempty(I), continue; end

            % 严格大圆距离筛 ≤ sth
            dkm = great_circle_km(latr(I), lonr(I), latr(J), lonr(J), prm.R);
            keep = dkm <= prm.sth + 1e-12;
            if ~any(keep), continue; end
            I = I(keep); J = J(keep); dkm = dkm(keep);

            lam = 0.5*(y(I) - y(J)).^2;

            tau = t_cur - t_prev;
            Dbuf{end+1} = as_col(dkm);
            Tbuf{end+1} = repmat(tau, numel(lam), 1);
            Lbuf{end+1} = as_col(lam);
            Pbuf{end+1} = repmat(t_prev, numel(lam), 1);
            Cbuf{end+1} = repmat(t_cur,  numel(lam), 1);
        end
    end

    % ---------- 输出拼装（安全） ----------
    if isempty(Dbuf) || isempty(Lbuf)
        S = [];
        return;
    end

    D = vcat_cols(Dbuf);%spatial distance
    T = vcat_cols(Tbuf);%time
    G = vcat_cols(Lbuf);%variogram value,0.5*xxx
    P = vcat_cols(Pbuf);%previous epoch time
    C = vcat_cols(Cbuf);%current epoch time

    % 进一步保险：长度一致性检查
    n = min([numel(D), numel(T), numel(G), numel(P), numel(C)]);
    if n==0, S=[]; return; end
    D=D(1:n); T=T(1:n); G=G(1:n); P=P(1:n); C=C(1:n);

    S = [D, T, G, P, C];
    if prm.use_single, S = single(S); end
end

% ====== 小工具：列向量化 ======
function v = as_col(v), v = reshape(v, [], 1); end

% ====== 小工具：安全拼接列向量 cell，自动剔空，确保列向量 ======
function out = vcat_cols(C)
    if isempty(C), out = []; return; end
    % 剔除空单元
    C = C(~cellfun('isempty', C));
    if isempty(C), out = []; return; end
    % 全部列化
    for i=1:numel(C), C{i} = reshape(C{i}, [], 1); end
    out = vertcat(C{:});
end

% ===== 工具：向量化大圆距离（km） =====
function dkm = great_circle_km(lat1,lon1,lat2,lon2,R)
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    s = sin(0.5*dlat).^2 + cos(lat1).*cos(lat2).*sin(0.5*dlon).^2;
    ang = 2*atan2( sqrt(max(s,0)), sqrt(max(1-s,0)) );
    dkm = R * ang;
end

% ===== 无 Statistics Toolbox 的邻居搜索替代 =====
function neigh = local_rangesearch_fallback(U, r_chord)
    n = size(U,1);
    neigh = cell(n,1);
    blk = 8000;  % 可按内存调
    for i1 = 1:blk:n
        i2 = min(n, i1+blk-1);
        D = 1 - U(i1:i2,:)*U.';          % 1 - dot ∈ [0,2]
        D = sqrt(max(0, 2*D));           % chord 距离
        for ii = 1:(i2-i1+1)
            js = find(D(ii,:) <= r_chord + 1e-12);
            neigh{i1+ii-1} = js(:).';
        end
    end
end

% ===== 参数解析（与你现有风格一致） =====
function X = parseargs(X,varargin)
    remaining = nargin-1; count = 1;
    fields = fieldnames(X);
    modified = zeros(size(fields));
    while remaining>=2
        fieldname = varargin{count};
        fieldind = find(strcmp(fieldname,fields));
        if ~isempty(fieldind)
            oldvalue = getfield(X,fieldname); %#ok
            newvalue = varargin{count+1};
            if iscell(oldvalue)
                if ~iscellstr(oldvalue), error('All allowed values for "%s" must be strings',fieldname); end
                if ~ischar(newvalue),   error('New value for "%s" must be a string',fieldname); end
                if isempty(find(strcmp(oldvalue,newvalue),1))
                    error('"%s" is not allowed for field "%s"',newvalue,fieldname);
                end
            elseif ~isempty(oldvalue)
                if ~strcmp(class(oldvalue),class(newvalue))
                    error('Cannot change class of field "%s" from "%s" to "%s"',fieldname,class(oldvalue),class(newvalue));
                elseif numel(oldvalue)==1 && numel(newvalue)~=1
                    error('New value for "%s" must be a scalar',fieldname);
                end
            end
            X = setfield(X,fieldname,newvalue); %#ok
            modified(fieldind) = 1;
        else
            error(['Not a valid field name: ' fieldname]);
        end
        remaining = remaining - 2; count = count + 2;
    end
    if remaining~=0, error('Odd number of arguments supplied. Name-value pairs required'); end
    notmodified = find(~modified);
    for i=1:length(notmodified)
        fieldname = fields{notmodified(i)};
        oldvalue = getfield(X,fieldname); %#ok
        if iscell(oldvalue)
            if ~iscellstr(oldvalue), error('All allowed values for "%s" must be strings',fieldname); end
            if isempty(oldvalue),    error('Empty cell array not allowed for field "%s"',fieldname); end
            X = setfield(X,fieldname,oldvalue{1}); %#ok
        end
    end
end
