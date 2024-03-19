clear 
clc
% DH = [0 0 0 0;
%       0 0 0.34 0;
%       0 -pi/2 0 0;
%       0 pi/2 0.4 0;
%       0 pi/2 0 0;
%       0 -pi/2 0.4 0;
%       0 -pi/2 0 0;
%       0 pi/2 0.126 0];
DH = [0 0 0 0;   % shifted DH
      0 0 0 0;
      0 -pi/2 0 0;
      0 pi/2 0.4 0;
      0 pi/2 0 0;
      0 -pi/2 0.4 0;
      0 -pi/2 0 0;
      0 pi/2 0 0];
numjoints = size(DH,1)-1;
%%
syms m [1 numjoints] real
syms cx cy cz ixx ixy ixz iyy iyz izz [numjoints 1] real
syms g real
for i = 1:numjoints
    COM(i,:) = [cx(i),cy(i),cz(i)]';
end
for i = 1:numjoints
    I(:,:,i) = [ixx(i),ixy(i),ixz(i);
                ixy(i), iyy(i), iyz(i);
                ixz(i), iyz(i), izz(i)]';
end
for i = 1:7
    phi_(i,:) = [m(i), m(i)*COM(i,:), I(1,1,i),I(1,2,i),I(1,3,i),I(2,2,i),I(2,3,i),I(3,3,i)];
end
phi = [phi_(1,:),phi_(2,:),phi_(3,:),phi_(4,:),phi_(5,:),phi_(6,:),phi_(7,:)]';
%% matlab model check
iiwa = loadrobot("kukaIiwa14");
iiwa.DataFormat = "row";
iiwa.Gravity = [0,0, -9.8];
d_pose = iiwa.homeConfiguration;
%% Calculate torques first on random configruations
ra = -2;rb = 2; rN=7;
% count = 0;
Nn = 100;
for i  = 1:Nn
    q(i,:) = ra + (rb-ra).*rand(1,rN);
    qdot(i,:) = rand(1,rN);
    qddot(i,:) = rand(1,rN);
    % show(iiwa,q(i,:))
    tau(i,:) = gravityTorque(iiwa,q(i,:));
    % tau(i,:) = inverseDynamics(iiwa,q(i,:),qdot(i,:),qddot(i,:));
    K{i,:} = RegressorMatrix(DH,[0,q(i,:)],zeros(1,8),zeros(1,8));
    % K{i,:} = RegressorMatrix(DH,[0,q(i,:)],[0,qdot(i,:)],[0,qddot(i,:)]);
    % tau_reg(i,:) = RegressorMatrix(DH,[0,q(i,:)],zeros(1,8),zeros(1,8))*phi;
    % pause
    % tau_v(count,:) = tau(i,:)';
    % count = count+7;
end
tau_v = reshape(tau',[i*7,1]);

Y = cell2mat(K);
%%
zero_columns = all(Y == 0);
non_iden = find(zero_columns);
Y_iden = Y;
% Y_iden(:, non_iden) = [];
%%
[Q, R, P] = qr(Y_iden);
b = rank(Y_iden);
P_p = P(:,1:b);
P_d = P(:,b+1:end);
R_p = R(1:b,1:b);
R_d = R(1:b,b+1:end);
K_factor = inv(R_p)*R_d;
beta = [P_p' + K_factor*P_d'];
base_params = beta*phi;
base_params_ = vpa(base_params);
threshold = 1e-10;
mapSymType(base_params_, 'vpareal', @(x) piecewise(abs(x)<=threshold, 0, x))
 
% Y_b = Y_iden*pinv(beta);
Y_b = Y_iden*P_p;
solved_Y_b = inv(Y_b'*Y_b)*Y_b';
params = solved_Y_b*tau_v
%%
tau_ver = Y_b*params;
tau_ver_mat = reshape(tau_ver,[7,Nn]);
tau_ver_mat = tau_ver_mat';
%%
for i =1:7
    plot(tau(:,i));
    hold on
    plot(tau_ver_mat(:,i))
    hold off
    figure
end
