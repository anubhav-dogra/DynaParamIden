function [base_params, beta, Y_b]= base_parameters(DH,q,qdot,qddot)
    %% initializing all inertial parameters
    numjoints = size(DH,1)-1;
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

    %%  generatoring Regressor Matrix
    Nn = size(q,1);
    for i  = 1:Nn
        K{i,:} = RegressorMatrix(DH,[0,q(i,:)],[0,qdot(i,:)],[0,qddot(i,:)]);
    end
       Y = cell2mat(K);
    %% Segmenting base parameters
    zero_columns = all(Y == 0);
    non_iden = find(zero_columns);
    Y_iden = Y;
    % Y_iden(:, non_iden) = [];
    %%
    [Q, R, P] = qr(Y_iden);
    b = rank(Y_iden);
    % p = size(R,2);
    P_p = P(:,1:b);
    P_d = P(:,b+1:end);
    R_p = R(1:b,1:b);
    R_d = R(1:b,b+1:end);
    K_factor = inv(R_p)*R_d;
    beta = P_p' + K_factor*P_d';
    base_params__ = beta*phi;
    base_params_ = vpa(base_params__);
    threshold = 1e-10;
    base_params = mapSymType(base_params_, 'vpareal', @(x) piecewise(abs(x)<=threshold, 0, x));
     
    % Y_b = Y_iden*pinv(beta);
    Y_b = Y_iden*P_p;
    % phi_s = pinv(beta)*base_params;

    % inv_map = P*[ones(b), -K_factor;
    %             zeros(p-b,b), ones(p-b)];
    % phi_d = P_d'*phi;
    % stand_params = inv_map*[base_params;phi_d];
end