% Function for 2D bin grouping
function [s_centers, t_unique, grouped_exp_grid] = bin_grouping_2d(s_lags, t_lags, exp_values, num_s_bins)
    
s_lags   = s_lags(:);
t_lags   = t_lags(:);
exp_values = exp_values(:);
% ---------- 过滤非法值 ----------
finite = isfinite(s_lags) & isfinite(t_lags) & isfinite(exp_values);
s_lags = s_lags(finite);
t_lags = t_lags(finite);
exp_values = exp_values(finite);

if isempty(s_lags)
    s_centers = []; t_unique = []; grouped_exp_grid = [];
    return;
end

% Determine unique temporal lags
    [t_unique, ~, t_group_idx] = unique(t_lags);
    num_t_unique = length(t_unique);

    % Bin grouping for s_lags into num_s_bins bins
    min_s = min(s_lags);
    max_s = max(s_lags);
    s_bin_edges = linspace(min_s-1e-10, max_s+1e-10, num_s_bins + 1);
    [~, ~, s_bin_idx] = histcounts(s_lags, s_bin_edges);
    
    % 只取落入 bin 的样本
    valid = s_bin_idx > 0;
    if ~any(valid) %空间量均为空的话不处理，返回空
        s_centers = []; t_unique = []; grouped_exp_grid = [];
        return;
    end
    s_bin_idx    = s_bin_idx(valid);
    t_group_idx  = t_group_idx(valid);
    exp_values   = exp_values(valid);

    
    % Compute mean exp_values for each 2D bin, using NaN for empty bins
    grouped_exp_grid = accumarray([t_group_idx, s_bin_idx], exp_values, [num_t_unique,num_s_bins], @mean, NaN);
    
    % Compute spatial bin centers
    s_centers = (s_bin_edges(1:end-1) + s_bin_edges(2:end)) / 2;
end