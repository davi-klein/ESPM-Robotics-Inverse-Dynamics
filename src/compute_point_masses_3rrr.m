function masses = compute_point_masses_3rrr(robot, kin)
% COMPUTE_POINT_MASSES_3RRR Calculates point masses kinematics for 3-RRR.
%
% Determines the global position, velocity, and acceleration of the discrete 
% point masses equimomental to the 6 leg links and the moving platform.
%
% Args:
%   robot (struct): Structural and inertial parameters.
%   kin (struct): Kinematic states (B, C, P positions, velocities, accels).
%
% Returns:
%   masses (struct): Kinematic data containing:
%       - p_O (3x3x7 double): Global pos of point masses [coords, points, body].
%       - vp_O (2x3x7 double): Global vel of point masses.
%       - ap_O (2x3x7 double): Global accels of point masses.
%         * Indices 1..6 represent the legs, index 7 represents the platform.

    % Define the equilateral triangle vertices in the local frame
    q = [ 0,         sqrt(6)/2, -sqrt(6)/2;
          sqrt(2),  -sqrt(2)/2, -sqrt(2)/2;
          1,         1,          1        ];

    % Pre-allocate memory for speed (7 bodies, 3 masses each, X/Y/1 coords)
    masses.p_O = zeros(3, 3, 7);
    masses.vp_O = zeros(2, 3, 7);
    masses.ap_O = zeros(2, 3, 7);

    %=====================================================================%
    %                       6 LEG LINKS KINEMATICS                        %
    %=====================================================================%
    for j = 1:3
        % --- Lower Legs (Links 1, 3, 5) attached to Base A ---
        idx_L = 2*j - 1;
        A = robot.base(j, :)';
        B = kin.B(:, j); vB = kin.vB(:, j); aB = kin.aB(:, j);
        L_L = robot.L(idx_L);

        D_L = sqrt((1 / robot.legs(idx_L).m) * robot.legs(idx_L).E);
        p_A = [1, 0, robot.legs(idx_L).xcm_joint; 0, 1, 0; 0, 0, 1] * (D_L * q);
        
        R_L = (1/L_L) * [(B(1)-A(1)), -(B(2)-A(2));
                         (B(2)-A(2)),  (B(1)-A(1))];
        
        masses.p_O(:,:,idx_L) = [R_L, A; 0, 0, 1] * p_A;
        
        dR_L = (1/L_L) * [vB(1), -vB(2); vB(2), vB(1)];
        ddR_L = (1/L_L) * [aB(1), -aB(2); aB(2), aB(1)];
        
        masses.vp_O(:,:,idx_L) = dR_L * p_A(1:2, :);
        masses.ap_O(:,:,idx_L) = ddR_L * p_A(1:2, :);

        % --- Upper Legs (Links 2, 4, 6) attached to Joint B ---
        idx_U = 2*j;
        C = kin.C(:, j); vC = kin.vC(:, j); aC = kin.aC(:, j);
        L_U = robot.L(idx_U);

        D_U = sqrt((1 / robot.legs(idx_U).m) * robot.legs(idx_U).E);
        p_B = [1, 0, robot.legs(idx_U).xcm_joint; 0, 1, 0; 0, 0, 1] * (D_U * q);
        
        R_U = (1/L_U) * [(C(1)-B(1)), -(C(2)-B(2));
                         (C(2)-B(2)),  (C(1)-B(1))];
                         
        masses.p_O(:,:,idx_U) = [R_U, B; 0, 0, 1] * p_B;
        
        dR_U = (1/L_U) * [(vC(1)-vB(1)), -(vC(2)-vB(2));
                          (vC(2)-vB(2)),  (vC(1)-vB(1))];
        ddR_U = (1/L_U) * [(aC(1)-aB(1)), -(aC(2)-aB(2));
                           (aC(2)-aB(2)),  (aC(1)-aB(1))];
                           
        masses.vp_O(:,:,idx_U) = repmat(vB, 1, 3) + dR_U * p_B(1:2, :);
        masses.ap_O(:,:,idx_U) = repmat(aB, 1, 3) + ddR_U * p_B(1:2, :);
    end

    %=====================================================================%
    %                  MOVING PLATFORM KINEMATICS (BODY 7)                %
    %=====================================================================%
    D_P = sqrt((1 / robot.platform.m) * robot.platform.E);
    pef = D_P * q; % Points in center of mass frame
    
    C3 = kin.C(:, 3);
    vC3 = kin.vC(:, 3);
    aC3 = kin.aC(:, 3);
    
    P = kin.P;
    vP = kin.vP;
    aP = kin.aP;
    
    % Kinematics matching original formulation based on distance C3 to P
    dpc3 = norm(C3 - P);
    RP_O = (1/dpc3) * [(C3(1)-P(1)), -(C3(2)-P(2));
                       (C3(2)-P(2)),  (C3(1)-P(1))];
                       
    masses.p_O(:,:,7) = [RP_O, P; 0, 0, 1] * pef;
    
    dRP_O = (1/dpc3) * [(vC3(1)-vP(1)), -(vC3(2)-vP(2));
                        (vC3(2)-vP(2)),  (vC3(1)-vP(1))];
    masses.vp_O(:,:,7) = repmat(vP, 1, 3) + dRP_O * pef(1:2, :);
    
    ddRP_O = (1/dpc3) * [(aC3(1)-aP(1)), -(aC3(2)-aP(2));
                         (aC3(2)-aP(2)),  (aC3(1)-aP(1))];
    masses.ap_O(:,:,7) = repmat(aP, 1, 3) + ddRP_O * pef(1:2, :);
end