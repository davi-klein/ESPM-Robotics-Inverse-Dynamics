function torques = compute_inverse_dynamics_2r(robot, kin, masses)
% COMPUTE_INVERSE_DYNAMICS Calculates required joint torques.
%
% Applies D'Alembert's principle with equimomental systems of point masses
% and graph theory (Davies' method) to compute required actuation torques
% for the manipulator.
%
% Args:
%   robot (struct): Structural and inertial parameters.
%   kin (struct): Kinematic states (B, C positions, etc.).
%   masses (struct): Kinematic states of equimomental point masses.
%
% Returns:
%   torques (struct): Computed torques and forces, containing:
%       - joint_A (double): Required torque at joint A (N.m).
%       - joint_B (double): Required torque at joint B (N.m).

    % 1. Inertial Moment and Forces Screws (D'Alembert)
    S = [0, -1; 1, 0];
    
    S_in1 = zeros(3, 1);
    S_in2 = zeros(3, 1);
    
    % Elegant summation over the 3 point masses for link 1
    for j = 1:3
        p1_j = masses.p1_O(1:2, j);
        ap1_j = masses.ap1_O(1:2, j);
        S_in1 = S_in1 + [(S * p1_j)' * ap1_j; ap1_j];
    end
    S_in1 = -(robot.m1 / 3) * S_in1;
    
    % Elegant summation over the 3 point masses for link 2
    for j = 1:3
        p2_j = masses.p2_O(1:2, j);
        ap2_j = masses.ap2_O(1:2, j);
        S_in2 = S_in2 + [(S * p2_j)' * ap2_j; ap2_j];
    end
    S_in2 = -(robot.m2 / 3) * S_in2;

    % 2. External Forces and Torques Screw
    tep = 0; fex = 0; fey = 0; % Assuming no environmental contact forces
    te = tep - fex * kin.C(2) + fey * kin.C(1);
    S_ex = [te; fex; fey];

    % 3. Weight Forces Screws (Requires CoM in inertial frame)
    x_a = 0; y_a = 0;
    x_b = kin.B(1); y_b = kin.B(2);
    x_c = kin.C(1); y_c = kin.C(2);
    
    RA_O = (1 / robot.L1) * [(x_b - x_a), -(y_b - y_a);
                             (y_b - y_a),  (x_b - x_a)];
    cm1_O = [RA_O, [x_a; y_a]; 0, 0, 1] * [robot.L1 / 2; 0; 1];
    S_W1 = [cm1_O(1) * robot.m1 * robot.g; 0; robot.m1 * robot.g];
    
    RB_O = (1 / robot.L2) * [(x_c - x_b), -(y_c - y_b);
                             (y_c - y_b),  (x_c - x_b)];
    cm2_O = [RB_O, kin.B; 0, 0, 1] * [robot.L2 / 2; 0; 1];
    S_W2 = [cm2_O(1) * robot.m2 * robot.g; 0; robot.m2 * robot.g];

    % 4. Primary Variables Unit Screws (Joint Constraints)
    s_AFx = [-y_a; 1; 0];
    s_AFy = [x_a; 0; 1];
    s_Atz = [1; 0; 0];
    
    s_BFx = [-y_b; 1; 0];
    s_BFy = [x_b; 0; 1];
    s_Btz = [1; 0; 0];

    % 5. Cutset Matrix (Graph Theory for Davies Method)
    % Columns: s_AFx, s_AFy, s_Atz, s_BFx, s_BFy, s_Btz, S_ex, S_in1, S_in2, S_W1, S_W2
    Q = [ 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1;   % Fundamental cut A
          0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 ]; % Fundamental cut B 

    % 6. Network Actions Matrix
    A_d = [s_AFx, s_AFy, s_Atz, s_BFx, s_BFy, s_Btz, S_ex, S_in1, S_in2, S_W1, S_W2];
    
    A_n = [A_d * diag(Q(1, :));
           A_d * diag(Q(2, :))];
           
    A_n_s = A_n(:, 1:6);  % Secondary variables (Unknowns)
    A_n_p = A_n(:, 7:11); % Primary variables (Knowns)
    psi_p = [1; 1; 1; 1; 1]; % Magnitudes of primary actions

    % 7. Solve for Unknown Torques/Forces
    Psi = A_n_s \ (-A_n_p * psi_p);
    
    % Extract the required joint torques
    torques.joint_A = Psi(3);
    torques.joint_B = Psi(6);
end