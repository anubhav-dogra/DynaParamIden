clc
clear
DH = [0 0 0.34 0;
      0 -90 0 0;
      0 90 0.4 0;
      0 90 0 0;
      0 -90 0.4 0;
      0 -90 0 0;
      0 90 0.126 0]; % coz the loop will start from i = 2
% DH = [0 0 0 0;
%       0 -pi/2 0 0;
%       0 pi/2 0.4 0;
%       0 pi/2 0 0;
%       0 -pi/2 0.4 0;
%       0 -pi/2 0 0;
%       0 pi/2 0 0];
numjoints = size(DH,1);
syms q qdot qddot m [1 numjoints] real
syms cx cy cz ixx ixy ixz iyy iyz izz [numjoints 1] real
syms g real
for i = 1:numjoints
    COM(:,i) = [cx(i),cy(i),cz(i)]';
end


%%
w = sym(zeros(3,numjoints+1));wdot = sym(zeros(3,numjoints+1));
vdot(:,1) = [0;0;-g];
z(:,1) = [0;0;1];
count = 1;
T1 = eye(4);
for i = 2:numjoints+1    % cant do i-1 as 0, therefore i-1 is 1 and i =  2. 
    T(:,:,i) = TM(DH(i-1,1),DH(i-1,2),DH(i-1,3),q(i-1)); % T(:,:,2) = is actually 0_T_1
    T1 = T1*T(:,:,i);
    frame_wrt_base(:,:,i-1) = T1;
    Z_i = [0 0 1]';
    % Z_i = frame_wrt_base(1:3,3,i-1);

    w(:,i) = T(1:3,1:3,i)'*w(:,i-1) + qdot(i-1)*Z_i;
    wdot(:,i) = T(1:3,1:3,i)'*wdot(:,i-1) + cross((T(1:3,1:3,i)'*w(:,i-1)),(qdot(i-1)*Z_i))+qddot(i-1)*Z_i;
    vdot(:,i) = T(1:3,1:3,i)'*(cross(wdot(:,i-1),T(1:3,4,i)) + cross(w(:,i-1),(cross(w(:,i-1),T(1:3,4,i))))+vdot(:,i-1));
    vcdot(:,i-1) = cross(wdot(:,i),COM(:,i-1)) + cross(w(:,i),(cross(w(:,i),COM(:,i-1)))) + vdot(:,i);
    COM_base(:,i-1) = frame_wrt_base(:,:,i-1)*[COM(:,i-1);1];
    U_l(i-1) = m(i-1)*[0 0 -g]*COM_base(1:3,i-1);
end
U = -sum(U_l);
% simplify(vpa(U))
G = gradient(U,q);
% G = vpa(G);
solvedG = subs(G,q,[0.0602    1.7022    0.7651   -1.6691   -0.6477   -1.8656    0.0078]);
solvedG = subs(solvedG,g,9.8);
solvedG = subs(solvedG,m,[4,4,3,2.7,1.7,1.8, 0.3]);
solvedG = subs(solvedG,cx,[0,-0.0003,0,0,-0.0001,0,0]');
solvedG = subs(solvedG,cy,[-0.0013,-0.0059,0.03,0.067,-0.021,-0.0006,0]');
solvedG = subs(solvedG,cz,[-0.0825,0.042,-0.0755,0.034,-0.1395,0.0004,0.02]');
%%

% subs(G,[q, g],[0.0602    1.7022    0.7651   -1.6691   -0.6477   -1.8656    0.0078, 9.8])
% tau = [0.0000   41.3264   14.1277    0.3259   -0.3426   -0.0000         0]
%% Transformation Matrix
function T = TM(a, alph, d, q)
T = [cos(q), -sin(q), 0, a;
        sin(q)*cosd(alph), cos(q)*cosd(alph), -sind(alph), -sind(alph)*d;
        sin(q)*sind(alph), cos(q)*sind(alph), cosd(alph), cosd(alph)*d;
        0,0,0,1];
end