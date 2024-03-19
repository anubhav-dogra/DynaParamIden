
tau_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_effort.txt");
tau_ext = readtable("/home/terabotics/work/DynaParamIden/data/external_torques.txt");
tau_t = table2array(tau_t);
tau_v = reshape(tau_t',[size(tau_t,1)*7,1]);
tau_ext = table2array(tau_ext);
tau_ext_v = reshape(tau_ext',[size(tau_ext,1)*7,1]);

q_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_positions.txt");
q = table2array(q_t);

qdot_t = readtable("/home/terabotics/work/DynaParamIden/data/joint_velocity.txt");
qdot = 0*table2array(qdot_t);
Nn = size(q,1);
qddot = zeros(Nn,7);
%%
DH = [0 0 0 0;   % shifted DH
      0 0 0 0;
      0 -pi/2 0 0;
      0 pi/2 0.4 0;
      0 pi/2 0 0;
      0 -pi/2 0.4 0;
      0 -pi/2 0 0;
      0 pi/2 0 0];

[base_params, Y_b]= base_parameters(DH,q,qdot,qddot);

solved_Y_b = inv(Y_b'*Y_b)*Y_b';
params = solved_Y_b*tau_v;
% params = solved_Y_b*(tau_v-tau_ext_v);
%%
tau_ver = Y_b*params;
tau_ver_mat = reshape(tau_ver,[7,Nn]);
tau_ver_mat = tau_ver_mat';

%%
for i =1:7
    plot(tau_t(:,i));
    hold on
    plot(tau_ver_mat(:,i))
    hold off
    figure
end

