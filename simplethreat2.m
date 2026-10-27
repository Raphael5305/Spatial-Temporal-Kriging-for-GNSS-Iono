function results=simplethreat2(R,N,IPPVD,samp_int,Rin_int,IGPs,Broadcast_GIVD, Broadcast_IGPs,dq)
% we need the S_PF model data,
if nargin < 9, dq = []; end

R_min=R(1);R_max=R(2);% for Hong Kong: 50 and 400
N_target=N(2);N_min=N(1); % in WAAS: 10 and 30
wgs84 = wgs84Ellipsoid('kilometer');
height=(350000)/1000;%unit:km
R_earth = 6378.137;  % 地球半径（单位：km）

PF=IPPVD; % in this place, IPPVD should be PF detredending or nontrending,SPF data

% 规范 IGPs 经度 [-180,180]
IGP_lat = IGPs(:,1);
IGP_lon = mod(IGPs(:,2) + 180, 360) - 180;

% 规范 Broadcast_IGPs 经度
B_lat = Broadcast_IGPs(:,1);
B_lon = mod(Broadcast_IGPs(:,2) + 180, 360) - 180;
Broadcast_IGPs = [B_lat, B_lon];

% 历元集合（按 samp_int 下采样）
all_epochs = unique(IPPVD(:,1), 'stable');
est_epochs = all_epochs(1:samp_int:numel(all_epochs));
num_epochs = numel(est_epochs);
num_igps   = numel(IGP_lat);

results = zeros(((R(2)-R(1))/Rin_int+1)* num_epochs*num_igps, 14);
% --- IGP → Broadcast_IGPs 映射（精确匹配；你需要最近邻可自行加上） ---
tol = 1e-8;
B_map = nan(num_igps,1);
for ig = 1:num_igps
    mask_lat = abs(B_lat - IGP_lat(ig)) <= tol;
    mask_lon = abs(B_lon - IGP_lon(ig)) <= tol;
    idx = find(mask_lat & mask_lon, 1, 'first');
    if ~isempty(idx)
        B_map(ig) = idx;
    end
end

% === Default UnK 模型（用于 UnK_threat 或 ST 缺省场景的排序距离核）===
defaultUnK.func        = @(b,h) b(2) * (1 - exp(-h ./ b(1)));
defaultUnK.sill        = 0.85;%1;           % total^2
defaultUnK.model       = 'exponential';
defaultUnK.nugget      = 0.05;%0.3 * 0.3;   % nom^2
defaultUnK.partialsill = 0.8;
defaultUnK.range       = 8000;        % km （空间相关距离，用于排序时的尺度）
defaultfit             = [defaultUnK.range, defaultUnK.partialsill, defaultUnK.nugget];

% PF 9 row，predict,
% chi,chi_irreg,R_irreg,sigma_inflated(IPP),sigma_undersampled,RCM, Rfit,
save_row=0;

for ie=1:num_epochs
    est_time = est_epochs(ie);
    time_idx = find(IPPVD(:,1) == est_time);            % 当前历元索引

    if isempty(time_idx)
        %若窗口空，直接跳过这个历元；若当前无历元也跳过
        continue;
    end

    % 日内 30s 索引 j（列1=00:00:00）
    jj= epoch_to_30s_index(est_time)-1;%判断哪个时间是对的


    truth_val2= NaN;
    if ~isnan(jj)
        m = B_map(1);
        if ~isnan(m), truth_val2 = Broadcast_GIVD(m, jj);end
    end

    % Dist=[];East2=[];North2=[];Up2=[];
    IPPs=PF(time_idx,3:4);
    lat1 = IGP_lat;  % IGP纬度
    lon1 = IGP_lon;  % IGP经度-180~180
    lat2 = IPPs(:,1);      % IPP纬度
    lon2 = IPPs(:,2);      % IPP经度

    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    Dist_full = R_earth * c;  % 所有IPP与IGP的距离

    [Dist,id]=sort(Dist_full);
    % num_IPPs = sum(Dist<=R_min);
    sort_num=time_idx;
    id=sort_num(id); %转化为在PF中的绝对索引

    selected = Dist <= R_max;
    idsearch_initial = id(selected); % 绝对索引
    % Dist_selected = Dist(selected);  % 对应的距离
    if numel(idsearch_initial)<N_min
        continue;
    end

    R_in_values = 0:Rin_int:2000;
    for r = 1:length(R_in_values)
        R_in = R_in_values(r);
        for Rw=(0.1:0.1:0.5).*2000
            R_out = R_in+Rw;
            PF_samples = PF(idsearch_initial, :);
            target_point = [est_time,1,lat1,lon1,truth_val2];
            result_once = depriveAndEstimate2(R,N,R_in, R_out, PF_samples, target_point, defaultUnK, defaultfit, wgs84, height, R_earth);
            if isnan(result_once)
                % 处理样本不足无法估计的情况
                continue;
            end
            save_row=save_row+1;
            % 存储或处理Zhat结果
            results(save_row,:)=[target_point result_once];
        end
    end

end
results = results(~all(results == 0, 2), :);%去除多余的行
end

%parsave([workspace '\GIVD dense\GIVD_' num2str(dense) '.mat'], 'PF_est', PF_est, 'OK_g_est', OK_g_est, 'UnK_g_est', UnK_g_est, 'AOK2_e_est', OK_, 'AOK2_g_est', AOK2_g_est, 'IDW_2_est', IDW_2_est);


% === 辅助函数：列1=00:00:00 → j = floor(sod/30)+1 ===
function j = epoch_to_30s_index(epoch_seconds)
if ~isfinite(epoch_seconds), j = NaN; return; end
sod = mod(epoch_seconds, 86400);
if sod < 0, sod = sod + 86400; end
j = floor(sod / 30) + 1; % 0→1, 30→2, ..., 86370→2880
end



