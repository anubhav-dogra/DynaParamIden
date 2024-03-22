clear
clc
% close all
tau_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_effort.txt");
tau_ext = readtable("/home/terabotics/work/DynaParamIden/data/external_torques.txt");
tau_cmd = readtable("/home/terabotics/work/DynaParamIden/data/tau_commanded.txt");
tau_t = table2array(tau_t);
tau_v = reshape(tau_t',[size(tau_t,1)*7,1]);
tau_cmd = table2array(tau_cmd);
tau_cmd_v = reshape(tau_cmd',[size(tau_cmd,1)*7,1]);
tau_ext = table2array(tau_ext);
tau_ext_v = reshape(tau_ext',[size(tau_ext,1)*7,1]);
tau_g  = tau_t-tau_ext;
tau_g_v = reshape(tau_g',[size(tau_g,1)*7,1]);
q_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_positions.txt");
q_cmd = readtable("/home/terabotics/work/DynaParamIden/data/q_commanded.txt");
q = table2array(q_t);
% q = q(1:10,:)
% q = table2array(q_cmd);
qdot_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_velocity1.txt");
qdot = 0*table2array(qdot_t);
%%
Nn = size(q,1);
qddot = zeros(Nn,7);
%%
DH = [0 0 0 0;   % shifted DH
      0 0 0 0;
      0 -pi/2 0 0;
      0 pi/2 0.42 0;
      0 pi/2 0 0;
      0 -pi/2 0.4 0;
      0 -pi/2 0 0;
      0 pi/2 0 0];

[base_params, beta, Y_b]= base_parameters(DH,q,qdot,qddot);

solved_Y_b = inv(Y_b'*Y_b)*Y_b';
% params = solved_Y_b*tau_v;
params = solved_Y_b*tau_g_v;
% params = solved_Y_b*(tau_cmd_v-tau_ext_v);

% params = solved_Y_b*(tau_v-tau_ext_v);
%%
tau_ver = Y_b*params;
tau_ver_mat = reshape(tau_ver,[7,Nn]);
tau_ver_mat = tau_ver_mat';

%%
figure
for i =1:7
    subplot(2,4,i)
    plot(tau_g(:,i),'LineWidth',2.0);
    hold on
    plot(tau_ver_mat(:,i),'LineWidth',2.0)
    title("Joint torque:" +num2str(i))
    set(gca,'LineWidth',0.75,'FontSize',18,'XMinorTick','on','YMinorTick','on','TickLength',[.01 0.1], 'XMinorGrid','on','YMinorGrid','on');
    hold off
    grid on
    box on
    % legend("Observed", "Estimated")
    % figure
end

%%
% iiwa = loadrobot("kukaIiwa14");
% iiwa.DataFormat = "row";
% iiwa.Gravity = [0,0, -9.8];
% d_pose = iiwa.homeConfiguration;
% % iiwa.Bodies{1,10}.Mass = 1
% count = 1
% for i  = 1:Nn
%     % show(iiwa,q(count,:))
%     tau_sim(i,:) = gravityTorque(iiwa,q(i,:));
% end
