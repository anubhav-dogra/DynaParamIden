clear 
clc
close all
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
      0 pi/2 0.42 0;
      0 pi/2 0 0;
      0 -pi/2 0.4 0;
      0 -pi/2 0 0;
      0 pi/2 0 0];
%% matlab model check
iiwa = loadrobot("kukaIiwa14");
iiwa.DataFormat = "row";
iiwa.Gravity = [0,0, -9.8];
d_pose = iiwa.homeConfiguration;
% iiwa.Bodies{1,10}.Mass = 1
%% Calculate torques first on random configruations
ra = -2;rb = 2; rN=7;
% count = 0;
Nn = 20;
rng(1,"twister");
for i  = 1:Nn
    q(i,:) = ra + (rb-ra).*rand(1,rN);
    qdot(i,:) = 0*rand(1,rN);
    qddot(i,:) = 0*rand(1,rN);
    % show(iiwa,q(i,:))
    hold on
    tau(i,:) = gravityTorque(iiwa,q(i,:));
end
tau_v = reshape(tau',[i*7,1]);
% tau_v = tau_v- (0.5 + (0.5+0.5).*rand(i*7,1));
[base_params, beta, Y_b]= base_parameters(DH,q,qdot,qddot);
disp(base_params)
solved_Y_b = inv(Y_b'*Y_b)*Y_b';
params = solved_Y_b*tau_v
%%
tau_ver = Y_b*params;
tau_ver_mat = reshape(tau_ver,[7,Nn]);
tau_ver_mat = tau_ver_mat';
delta_m7 = (abs(params(7))-1.4956)/0.4
delta_m7_ = (abs(params(9))-3.9650)/0.4
%%
for i =1:7
    plot(tau(:,i));
    hold on
    plot(tau_ver_mat(:,i))
    hold off
    figure
end
