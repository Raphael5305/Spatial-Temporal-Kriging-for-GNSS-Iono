% Helper function for bin grouping
function [bins, centers, grouped_values] = bin_grouping(lags, values, num_bins)
   % 仅返回“非空 bin”的中心与均值，用于后续 PSO 拟合
% 输出均为列向量；不包含任何 NaN

    % --- 基本检查 ---
    if nargin < 3 || ~isscalar(num_bins) || num_bins < 1
        error('num_bins must be a positive scalar.');
    end
    if isempty(lags) || isempty(values)
        bins=[];centers = []; grouped_values = [];
        return;
    end
    if numel(lags) ~= numel(values)
        error('lags and values must have the same length.');
    end

    % --- 过滤非法值 ---
    mask = isfinite(lags) & isfinite(values);
    lags = lags(mask); values = values(mask);
    if isempty(lags)
        bins=[];centers = []; grouped_values = [];
        return;
    end

 % --- 退化情形：所有 lags 相等 ---
    min_lag = min(lags); max_lag = max(lags);
    if min_lag == max_lag
        bins=min_lag(:);
        centers = min_lag(:);              % 1x1 → 列
        grouped_values = mean(values(:));  % 1x1 → 标量
        grouped_values = grouped_values(:);
        return;
    end
%正常分组
    bin_edges = linspace(min_lag-1e-10, max_lag+1e-10, num_bins+1);
    [~, ~, bin_idx] = histcounts(lags, bin_edges);
    
    % Compute bin centers
    centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    % Compute mean values per bin and count assignments
    valid_idx = bin_idx > 0; % Only consider valid bin assignments
    grouped_values = accumarray(bin_idx(valid_idx), values(valid_idx), [num_bins, 1], @mean, NaN);
    
    % Identify non-empty bins (where grouped_values is not NaN)
    non_empty = ~isnan(grouped_values);
    
    % Filter outputs to include only non-empty bins
    grouped_values = grouped_values(non_empty);
    centers = centers(non_empty)';
    
    % Adjust bin_edges to include only edges for non-empty bins
    % Include the left edge of the first non-empty bin and the right edge of the last
    non_empty_indices = find(non_empty);
    if ~isempty(non_empty_indices)
        bins = bin_edges([non_empty_indices; non_empty_indices(end)+1]);
    else
        bins = [];
    end
end