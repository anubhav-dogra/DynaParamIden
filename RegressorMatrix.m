function K = RegressorMatrix(DH,q,qdot,qddot)
    numjoints = size(DH,1)-1;
    w(:,1) = [0;0;0];wdot(:,1) = [0;0;0];
    vdot(:,1) = [0;0;-9.80];
    z(:,1) = [0;0;1];
    Xeye = eye(6);
    count = 1;
    T1 = eye(4);
    for i = 2:numjoints+1    % cant do i-1 as 0, therefore i-1 is 1 and i =  2. 
        T(:,:,i) = TM(DH(i,1),DH(i,2),DH(i,3),q(i)); % T(:,:,2) = is actually 0_T_1
        T1 = T1*T(:,:,i);
        frame_wrt_base(:,:,count) = T1;
        Z_i = [0 0 1]';
        % Z_i = frame_wrt_base(1:3,3,i-1);
        w(:,i) = T(1:3,1:3,i)'*w(:,i-1) + qdot(i)*Z_i;
        wdot(:,i) = T(1:3,1:3,i)'*wdot(:,i-1) + cross((T(1:3,1:3,i)'*w(:,i-1)),(qdot(i)*Z_i))+qddot(i)*Z_i;
        vdot(:,i) = T(1:3,1:3,i)'*(cross(wdot(:,i-1),T(1:3,4,i)) + cross(w(:,i-1),(cross(w(:,i-1),T(1:3,4,i))))+vdot(:,i-1));
        A_n(:,:,count) = getAnMatrix(vdot(:,i),wdot(:,i),w(:,i));
        count = count+1;
    end
    
    for i = 1:numjoints-1
        X(:,:,i) = getX(T(:,:,i+2)); % coz 1_X_2 = 1_T_2
    end
    
    for i = 1:numjoints
        X_temp =Xeye;
        Z_i = [0 0 1];
        % Z_i = frame_wrt_base(1:3,3,i)';
        for j = 1:numjoints
            if i == j
                A{i,j} = eye(6)*A_n(:,:,i);
                K_{i,j} = [Z_i,0,0,0]*A{i,j};
                continue
    
            elseif i > j
                A{i,j} = zeros(6,10);
                K_{i,j} = [Z_i,0,0,0]*A{i,j};
                continue
    
            else
                X_temp= X_temp*X(:,:,j-1);
                
            end
            A{i,j} = X_temp*A_n(:,:,j);
            K_{i,j} = [Z_i,0,0,0]*A{i,j};
    
        end
        
    end
    K = cell2mat(K_);

% for i = 1:numjoints
%     phi_(i,:) = [M(i), M(i)*COM(i,1:3), I_l(1,1,i),I_l(1,2,i),I_l(1,3,i),I_l(2,2,i),I_l(2,3,i),I_l(3,3,i)];
% end
% phi = [phi_(1,:),phi_(2,:),phi_(3,:),phi_(4,:),phi_(5,:),phi_(6,:),phi_(7,:)]';
% % phi = [phi_(1,:),phi_(2,:),phi_(3,:)]';
% tau = K*phi;
end
% disp(tau)
%% spatial force transform
function X = getX(T)
    a_R_b = T(1:3,1:3);
    P = T(1:3,4);
    % X = [a_R_b, zeros(3);
    %     getskew(P)*a_R_b, a_R_b]; % for [force;torques] format
    X = [a_R_b -getskew(P)*a_R_b,;
        zeros(3) a_R_b];% for [torques;forces] format
end

%% spatial force matrix
function A_n_Matrix = getAnMatrix(d0ddot,wdot, w)
    A_n_Matrix = [zeros(3,1), -getskew(d0ddot), getL(wdot)+getskew(w)*getL(w);
                d0ddot, getskew(wdot)+getskew(w)*getskew(w) zeros(3,6)];
end


%% skew matrix
function skew_matrix = getskew(vector)
    skew_matrix = [0, -vector(3), vector(2);
                   vector(3), 0, -vector(1);
                   -vector(2), vector(1), 0];
end


%% L matrix
function L = getL(w)
    L = [w(1), w(2), w(3),0, 0, 0;
         0, w(1), 0, w(2), w(3), 0;
         0, 0, w(1), 0, w(2), w(3)];
end

%% Transformation Matrix
function T = TM(a, alph, d, q)
T = [cos(q), -sin(q), 0, a;
        sin(q)*cos(alph), cos(q)*cos(alph), -sin(alph), -sin(alph)*d;
        sin(q)*sin(alph), cos(q)*sin(alph), cos(alph), cos(alph)*d;
        0,0,0,1];
end
