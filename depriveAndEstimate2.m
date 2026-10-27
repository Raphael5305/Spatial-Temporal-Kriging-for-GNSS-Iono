function Zhat = depriveAndEstimate2(R,N,Rin, Rout, PF_samples, target_point, vstruct_UnK, UnK_fit, wgs84, height, Rearth)
    % depriveAndEstimate - 执行数据剔除和universal Kriging估计
    %
    % 输入:
    %   Rin           - 数据剔除的内半径 (km)
    %   Rout          - 数据剔除的外半径 (km)
    %   PF_samples    - 样本点矩阵 (IPP数据: 纬度, 经度, 值等)
    %   target_point  - 待估计的目标点 (IGP: 时间, 纬度, 经度等)
    %   vstruct_UnK   - universal Kriging的变差函数结构
    %   UnK_fit       - Kriging拟合参数
    %   N_min         - 估计所需的最小样本数
    %   wgs84         - WGS84椭球模型
    %   height        - 坐标转换的高度 (km)
    %   R             - 地球半径 (km)
    %
    % 输出:
    %   Zhat          - 目标点的估计值，若估计失败则返回NaN

    
    Rmin=R(1);
    Rmax=R(2);

    Nmin=N(1);Ntarget=N(2);
    % 提取目标点坐标
    lat1 = target_point(3); % IGP纬度
    lon1 = target_point(4); % IGP经度

    % 提取样本点坐标
    lat2 = PF_samples(:, 3); % IPP纬度
    lon2 = PF_samples(:, 4); % IPP经度

    % 使用Haversine公式计算距离
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    Dist = Rearth * c; % 距离 (km)

    % 创建掩码，剔除Rin < Dist <= Rout的样本
    mask = (Dist <= min(Rin,Rout))| (Dist > max(Rin,Rout));
    % 选择剩余样本
    PF_remaining = PF_samples(mask, :);
    Dist_remaining = Dist(mask);

    % 检查剩余样本数量
     if size(PF_remaining, 1) < Nmin
        Zhat = NaN; % 样本不足
        return;
    elseif size(PF_remaining, 1) > Ntarget
        PF_remaining = PF_remaining(1:Ntarget, :);
        Dist_remaining = Dist_remaining(1:Ntarget);
        Rfitt=max(max(Dist_remaining),Rmin);
     else % Nmin to Ntarget
         Rfitt=Rmax;
     end

    % 准备Kriging：转换为ENU坐标
    East2 = zeros(size(PF_remaining, 1), 1);
    North2 = zeros(size(PF_remaining, 1), 1);
    for k = 1:size(PF_remaining, 1)
        [East, North, ~] = geodetic2enu(...
            PF_remaining(k, 3), PF_remaining(k, 4), height, ...
            lat1, lon1, height, wgs84, 'degrees');
        East2(k) = East;
        North2(k) = North;
    end

    % 准备Kriging的距离矩阵
    dist_n0 = distance('gc', PF_remaining(:, 3), PF_remaining(:, 4), lat1, lon1, Rearth);

    num_samples = size(PF_remaining, 1);
    [M, N] = meshgrid(1:num_samples);
    lat_m = PF_remaining(M(:), 3);
    lon_m = PF_remaining(M(:), 4);
    lat_n = PF_remaining(N(:), 3);
    lon_n = PF_remaining(N(:), 4);
    dist_nn = reshape(distance('gc', lat_m, lon_m, lat_n, lon_n,Rearth), num_samples, num_samples);

    % RCM calculation,
    lati_mean=mean(PF_remaining(:, 3)); % average latitude and longitude
    long_mean=mean(PF_remaining(:, 4));
    Rcentroid1=distance('gc',lati_mean,long_mean,lat1,lon1,6378.137);
    RCM1=Rcentroid1/Rfitt;

    data_values = PF_remaining(:, 5); % 假设第5列为待估计值
    [X, Y] = meshgrid(0); % 目标点相对位置
    UnK_result=Unk_threat(vstruct_UnK, East2, North2, data_values, X, Y, UnK_fit, dist_nn, dist_n0);
    %zi,chi,chi_irreg,R_irreg,sigma_IPP,sigma_inflated
    a=(UnK_result(1)-target_point(5))^2/(5.33^2)-UnK_result(6);
    sigma_undersampled=sqrt(max(a,0));
    Zhat = [UnK_result Rfitt RCM1 sigma_undersampled];
end