function [results,condd] = STkriging(model,x,y,t_n0,z,xi,yi,ti,dis_nn,dis_n0,t_nn,obs)

%Changed Author: REN HUI repaired to fullfill the GPS measurements and correction 18th/June/2025
% RH Email: 22039387r@connect.polyu.hk
% z is the observations of sampling data

% size of input arguments
sizest = size(xi);
numobs = numel(x);

% force column vectors
xi = xi(:);
yi = yi(:);
ti = ti(:);
x  = x(:);
y  = y(:);
t  = t_n0(:);
z  = z(:);
%output weight


% now calculate the matrix with variogram values 
A = stpacked(model,dis_nn,t_nn);%含有nugget
C=-(A-model.sill);

% 根据条件数决定是否添加扰动项
[W,C_reg]=QRinv(C,2e5);
condd(1)=cond(C);
% 计算 W = inv(C_reg)，使用 QR 分解提高稳定性
%W=inv(C)
%% Check if t_n0 is all zeros or all identical
if all(t_n0 == 0) || all(t_n0 == t_n0(1))
    % Special case: remove time trend
    G = [ones(numobs, 1), x, y]; % n x 3
    num_trends = 3;
else
    % Normal case: include time trend
    G = [ones(numobs, 1), x, y, t_n0]; % n x 4
    num_trends = 4;
end
% Expand A matrix
A_exp = [[C_reg, G]; [G', zeros(num_trends, num_trends)]];

% A is often very badly conditioned. Hence we use the Pseudo-Inverse for
% solving the equations
[A_inv, ~] = QRinv(A_exp, 1e7);
condd(2)=cond(A_exp);
% we also need to expand z
z_exp = [z; zeros(num_trends, 1)];%RH: 4 constrants, 

% expand b with ones
b = [stpacked(model,dis_n0,t_n0)];
b = model.sill-b; % cpvariance matrix
b = [b; 1; zeros(num_trends-1, 1)];%changed by RH, as the extra constraints

%% solve system and output chi-square statistic
lambda = A_inv*b;% b鏄彉寮傚嚱鏁板?嘉?(hi0)
% estimate zi
lambda1=lambda(1:numobs);
z1 = lambda1'*z;
residual1=abs(z1-obs);
%% To avoid the lambda too big
%Method 2
    H = C_reg;
    f = -b(1:numobs);
    Aeq = G';
    beq = [1; zeros(num_trends-1, 1)];
    lb =-1* ones(numobs, 1); % 用户指定下界
    ub = 1* ones(numobs, 1); % 用户指定上界
    options = optimoptions('quadprog', 'MaxIterations', 50000,'Display', 'off');
    [lambda, fval, exitflag, output] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
    lambda2 = lambda / sum(lambda);  % 归一化
    z2 = lambda2' * z;  % 更新 zi
    residual2=abs(z2-obs);

if residual1<=residual2
    lambda=lambda1;
    zi=z1;
else 
    lambda=lambda2;
    zi=z2;
end
lambda3=ones(numobs,1)./numobs;
z3=lambda3' * z;
residual3=abs(z3-obs);
%% REN HUI
%calculate chi-square statistic,should be a constant
%matrix_C=% C noise matrix should be the constant part of varigriam function,size:N*N
%% Calculate Q and chi-square statistic
GWG = G' * W * G; % num_trends x num_trends
[GWG_inv, ~] = QRinv(GWG, 1e12); % Use QRinv for stability
condd(3)=cond(GWG);
Q = (W * G) * GWG_inv; % n x 4
chi=max(z(1:numobs)'*(W-Q*G'*W)*z(1:numobs),0);% chi-square
chi_number=numobs-num_trends;%random variable 
chi_irreg=chi/(chi2inv(0.999,chi_number));% default Pfa=10^-3;
R_irreg=max(1,chi2inv(0.999,chi_number)/chi2inv(0.001,chi_number));%平方项
 % lambda2= (W-Q*G'*W)*(b(1:numobs))+Q(:,1);%sparks weights, the same result

%% REN HUI: Kriging estimation error calculation
 sigma_IPP= max(0,lambda'*C*lambda-2*lambda'*b(1:numobs)+(model.sill));
% sigma_IPP=model.sill-lambda'*b; %blanch solution, the same result
sigma_inflated=R_irreg*sigma_IPP; %平方项

% zi = reshape(zi,sizest);
results=[zi,chi,chi_irreg,R_irreg,sigma_IPP,sigma_inflated];
% 检查非凸或未收敛
% if (residual3<=min([residual1,residual2]))
%     results=[];condd=[];
% end

