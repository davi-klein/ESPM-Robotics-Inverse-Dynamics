function torques = compute_inverse_dynamics_3rrr(robot, kin, masses)
% COMPUTE_INVERSE_DYNAMICS_3RRR Calculates required torques for the 3-RRR.
%
% Applies D'Alembert's principle via equimomental point masses and Davies' 
% method (graph theory) to compute the actuation torques for the 3 active joints.
%
% Args:
%   robot (struct): Structural and inertial parameters.
%   kin (struct): Kinematic states (B, C, P positions).
%   masses (struct): Kinematic states of the 7 equimomental systems.
%
% Returns:
%   torques (struct): Computed required torques:
%       - joint_A1, joint_A2, joint_A3 (double): Actuation torques (N.m).

    S = [0, -1; 1, 0];

    %=====================================================================%
    %               1. INERTIAL & WEIGHT SCREWS (7 BODIES)                %
    %=====================================================================%
    S_in = zeros(3, 7);
    S_W = zeros(3, 7);
    
    % Loop through the 6 legs
    for i = 1:6
        m_link = robot.legs(i).m;
        
        % D'Alembert Inertial Screws
        sum_in = zeros(3,1);
        for j = 1:3
            p_j = masses.p_O(1:2, j, i);
            ap_j = masses.ap_O(1:2, j, i);
            sum_in = sum_in + [(S * p_j)' * ap_j; ap_j];
        end
        S_in(:, i) = -(m_link / 3) * sum_in;
        
        % Weight Screws (Using equimomental CoM property: average of points)
        cm_O = mean(masses.p_O(1:2, :, i), 2);
        S_W(:, i) = [cm_O(1) * m_link * robot.g; 0; m_link * robot.g];
    end
    
    % End-Effector Platform (Body 7)
    m_plat = robot.platform.m;
    sum_in_plat = zeros(3,1);
    for j = 1:3
        p_j = masses.p_O(1:2, j, 7);
        ap_j = masses.ap_O(1:2, j, 7);
        sum_in_plat = sum_in_plat + [(S * p_j)' * ap_j; ap_j];
    end
    S_in(:, 7) = -(m_plat / 3) * sum_in_plat;
    
    cm_plat_O = mean(masses.p_O(1:2, :, 7), 2);
    S_W(:, 7) = [cm_plat_O(1) * m_plat * robot.g; 0; m_plat * robot.g];

    %=====================================================================%
    %               2. EXTERNAL FORCES & KINEMATIC SCREWS                 %
    %=====================================================================%
    tep = 0; fex = 0; fey = 0;
    te = tep - fex * kin.P(2) + fey * kin.P(1);
    S_ex = [te; fex; fey];

    % Programmatic generation of primary constraint screws
    s_A_Fx = zeros(3,3); s_A_Fy = zeros(3,3); s_A_tz = repmat([1;0;0], 1, 3);
    s_B_Fx = zeros(3,3); s_B_Fy = zeros(3,3);
    s_C_Fx = zeros(3,3); s_C_Fy = zeros(3,3);
    
    for k = 1:3
        % Joints A (Base)
        s_A_Fx(:, k) = [-robot.base(k, 2); 1; 0];
        s_A_Fy(:, k) = [ robot.base(k, 1); 0; 1];
        % Joints B (Active to Passive)
        s_B_Fx(:, k) = [-kin.B(2, k); 1; 0];
        s_B_Fy(:, k) = [ kin.B(1, k); 0; 1];
        % Joints C (Passive to Platform)
        s_C_Fx(:, k) = [-kin.C(2, k); 1; 0];
        s_C_Fy(:, k) = [ kin.C(1, k); 0; 1];
    end

    % Action Vectors Network (Ad) - Flattening the generated screws
    A_d = [s_A_Fx(:,1), s_A_Fy(:,1), s_A_tz(:,1), ...
           s_A_Fx(:,2), s_A_Fy(:,2), s_A_tz(:,2), ...
           s_A_Fx(:,3), s_A_Fy(:,3), s_A_tz(:,3), ...
           s_B_Fx(:,1), s_B_Fy(:,1), s_B_Fx(:,2), s_B_Fy(:,2), s_B_Fx(:,3), s_B_Fy(:,3), ...
           s_C_Fx(:,1), s_C_Fy(:,1), s_C_Fx(:,2), s_C_Fy(:,2), s_C_Fx(:,3), s_C_Fy(:,3), ...
           S_ex, S_in, S_W]; % S_in and S_W expand to their 7 columns each

    %=====================================================================%
    %               3. DAVIES METHOD (GRAPH CUTSET MATRIX)                %
    %=====================================================================%
    
    % Preserved exact topological cuts from original research
    Q = [ 1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
          1  1  1  0  0  0  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  1  1  1  0  1  1  1  1  1  1  0  1  1  1  1;
          1  1  1  0  0  0  1  1  1  1  1  0  0  0  0  0  0  1  1  0  0  1  1  1  0  0  1  1  1  1  1  0  0  1  1  1;
         -1 -1 -1  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0 -1  0  0  0  0  0  0;
         -1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0 -1 -1  0  0  0  0  0 -1 -1  0  0  0  0  0;
          0  0  0  0  0  0 -1 -1 -1  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0 -1  0  0;
          0  0  0  0  0  0 -1 -1 -1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0 -1 -1  0  0  0  0  0 -1 -1  0];

    % Network Actions Matrix
    A_n = [A_d * diag(Q(1,:));
           A_d * diag(Q(2,:));
           A_d * diag(Q(3,:));
           A_d * diag(Q(4,:));
           A_d * diag(Q(5,:));
           A_d * diag(Q(6,:));
           A_d * diag(Q(7,:))];

    % Secondary variables (first 21 columns) & Primary variables
    A_n_s = A_n(:, 1:21);
    A_n_p = A_n(:, 22:36);
    
    psi_p = ones(15, 1); % Magnitudes of primary actions (ex, in1..7, w1..7)

    % Linear solver for the unknown constraints and torques
    Psi = A_n_s \ (-A_n_p * psi_p);
    
    % The required torques correspond to the active joints z-axis (columns 3, 6, 9)
    torques.joint_A1 = Psi(3);
    torques.joint_A2 = Psi(6);
    torques.joint_A3 = Psi(9);
end