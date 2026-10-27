function S = variogram_temporal2(data, varargin)
    %计算不同 epoch 之间空间接近点对的时序变差函数值（暴力匹配版）
    % 输入：data = [lat, lon, epoch, y]
    % 输出：S 包含时间差、变差值和 epoch 对的表格

    % 设置默认参数
    params.sth = 5; % 距离阈值，单位：公里
    params.n = 10;    % 时间窗口，300秒
    params = parseargs(params, varargin{:});

    % 移除 NaN
    valid_idx = ~any(isnan(data), 2);
    data = data(valid_idx, :);

    % 获取并排序唯一 epoch
    epochs = unique(data(:, 3));
    epochs = sort(epochs);
    n_epochs = numel(epochs);

    % 地球半径（单位：公里）
    R = 6378.137;
    % data lack, only 1 epoch
    if n_epochs==1
        disp(['only 1 epoch at' num2str(epochs(1))]);
        S=[];
        return;
    end

   % 缓冲
    Tbuf = {}; Lbuf = {}; Pbuf = {}; Cbuf = {};

    for k = 2:n_epochs
        current_epoch = epochs(k);
        A = data(data(:,3)==current_epoch,:);
        if isempty(A), continue; end
        % 前150s窗口内的历元
        m_start = find(epochs >= current_epoch - 30*params.n, 1, 'first');

        for m = m_start:(k-1)
            prev_epoch = epochs(m);
            B = data(data(:,3)==prev_epoch,:);
            if isempty(B), continue; end

            for i = 1:size(A,1)
                % 用 distance 算大圆距离
                d = distance('gc', A(i,1), A(i,2), B(:,1), B(:,2), R);

                valid = d <= params.sth;
                if any(valid)
                    time_diff = current_epoch - prev_epoch;
                    lam = 0.5 * (A(i,4) - B(valid,4)).^2;

                    Tbuf{end+1} = repmat(time_diff, nnz(valid), 1);
                    Lbuf{end+1} = lam;
                    Pbuf{end+1} = repmat(prev_epoch, nnz(valid), 1);
                    Cbuf{end+1} = repmat(current_epoch, nnz(valid), 1);
                end
            end
        end
    end

    % 输出
    if isempty(Lbuf)
        S = [];
    else
        TemporalDistance = vertcat(Tbuf{:});
        Variogram        = vertcat(Lbuf{:});
        PreviousEpoch    = vertcat(Pbuf{:});
        CurrentEpoch     = vertcat(Cbuf{:});
        S = [TemporalDistance, Variogram, PreviousEpoch, CurrentEpoch];
    end
end

%% Subfunction: parseargs (modified to handle 'proximity_threshold')
function X = parseargs(X, varargin)
    remaining = nargin - 1; % Number of arguments other than X
    count = 1;
    fields = fieldnames(X);
    modified = zeros(size(fields));
    
    % Process input arguments two at a time
    while remaining >= 2
        fieldname = varargin{count};
        fieldind = find(strcmp(fieldname, fields));
        if ~isempty(fieldind)
            oldvalue = getfield(X, fieldname); %#ok
            newvalue = varargin{count + 1};
            % Handle 'proximity_threshold' specifically
            if strcmp(fieldname, 'proximity_threshold')
                if isnumeric(newvalue) && isscalar(newvalue) && newvalue > 0
                    X = setfield(X, fieldname, newvalue); %#ok
                    modified(fieldind) = 1;
                else
                    error('Threshold must be a positive scalar number');
                end
            else
                % Handle other fields (if any)
                if iscell(oldvalue)
                    if ~iscellstr(oldvalue)
                        error(sprintf('All allowed values for "%s" must be strings', fieldname)); %#ok
                    end
                    if ~ischar(newvalue)
                        error(sprintf('New value for "%s" must be a string', fieldname)); %#ok
                    end
                    if isempty(find(strcmp(oldvalue, newvalue))) %#ok
                        error(sprintf('"%s" is not allowed for field "%s"', newvalue, fieldname)); %#ok
                    end
                elseif ~isempty(oldvalue)
                    if ~strcmp(class(oldvalue), class(newvalue))
                        error(sprintf('Cannot change class of field "%s" from "%s" to "%s"', ...
                            fieldname, class(oldvalue), class(newvalue))); %#ok
                    elseif numel(oldvalue) == 1 && numel(newvalue) ~= 1
                        error(sprintf('New value for "%s" must be a scalar', fieldname)); %#ok
                    end
                end
                X = setfield(X, fieldname, newvalue); %#ok
                modified(fieldind) = 1;
            end
        else
            error(['Not a valid field name: ' fieldname]);
        end
        remaining = remaining - 2;
        count = count + 2;
    end
    
    % Check for odd number of arguments
    if remaining ~= 0
        error('Odd number of arguments supplied. Name-value pairs required');
    end
    
    % Handle unmodified cell arrays
    notmodified = find(~modified);
    for i = 1:length(notmodified)
        fieldname = fields{notmodified(i)};
        oldvalue = getfield(X, fieldname); %#ok
        if iscell(oldvalue)
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings', fieldname)); %#ok
            elseif isempty(oldvalue)
                error(sprintf('Empty cell array not allowed for field "%s"', fieldname)); %#ok
            end
            X = setfield(X, fieldname, oldvalue{1}); %#ok
        end
    end
end