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
%% matlab model check
iiwa = loadrobot("kukaIiwa14");
iiwa.DataFormat = "row";
iiwa.Gravity = [0,0, -9.8];
d_pose = iiwa.homeConfiguration;

%% Calculate torques first on random configruations
ra = -2;rb = 2; rN=7;
% count = 0;
Nn = 11;
for i  = 1:Nn
    q(i,:) = ra + (rb-ra).*rand(1,rN);
    qdot(i,:) = 0*rand(1,rN);
    qddot(i,:) = 0*rand(1,rN);
    % show(iiwa,q(i,:))
    tau(i,:) = gravityTorque(iiwa,q(i,:));
end
tau_v = reshape(tau',[i*7,1]);

[base_params, Y_b]= base_parameters(DH,q,qdot,qddot);

solved_Y_b = inv(Y_b'*Y_b)*Y_b';
params = solved_Y_b*tau_v;
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
