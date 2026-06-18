function masses = compute_point_masses(robot, kin)
% COMPUTE_POINT_MASSES Calculates the kinematics of equimomental point masses.
%
% This function determines the global position, velocity, and acceleration
% of the discrete point masses that are equimomental to the robot links.
% It uses non-rigid transformations and the kinematic states computed previously.
%
% Args:
%   robot (struct): The robot structural and inertial parameters.
%   kin (struct): Kinematic data of the robot joints (B, C) and their 
%                 velocities/accelerations, provided by compute_kinematics.
%
% Returns:
%   masses (struct): Kinematic data of the point masses including:
%       - p1_O, p2_O (2x3 double): Global positions of the point masses.
%       - vp1_O, vp2_O (2x3 double): Global velocities of the point masses.
%       - ap1_O, ap2_O (2x3 double): Global accelerations of the point masses.

    % Define base joint A (origin of the inertial reference frame)
    x_a = 0; 
    y_a = 0;
    
    % Unpack kinematic struct for cleaner code
    x_b = kin.B(1); y_b = kin.B(2);
    x_c = kin.C(1); y_c = kin.C(2);
    
    vx_b = kin.vB(1); vy_b = kin.vB(2);
    vx_c = kin.vC(1); vy_c = kin.vC(2);
    
    ax_b = kin.aB(1); ay_b = kin.aB(2);
    ax_c = kin.aC(1); ay_c = kin.aC(2);

    % Define the equilateral triangle vertices in the local frame
    q = [ 0,         sqrt(6)/2, -sqrt(6)/2;
          sqrt(2),  -sqrt(2)/2, -sqrt(2)/2;
          1,         1,          1        ];

    %=====================================================================%
    %                       LINK 1 POINT MASSES                           %
    %=====================================================================%
    
    % Non-rigid transformation to find points p1 in the CoM reference frame
    D1 = sqrt((1 / robot.m1) * robot.E1);
    p1 = D1 * q;
    
    % Center of mass coordinates of Link 1 relative to Joint A
    xcm1_A = robot.L1 / 2;
    ycm1_A = 0;
    
    % Homogeneous transformation from CoM to Joint A
    H1 = [1, 0, xcm1_A;
          0, 1, ycm1_A;
          0, 0, 1];
    p1_A = H1 * p1;
    
    % Transformation from Joint A to Inertial Frame (O)
    RA_O = (1 / robot.L1) * [(x_b - x_a), -(y_b - y_a);
                             (y_b - y_a),  (x_b - x_a)];
    HA_O = [RA_O, [x_a; y_a];
            0, 0, 1];
            
    masses.p1_O = HA_O * p1_A;
    
    % Velocity Analysis for Link 1 Point Masses
    dRA_O = (1 / robot.L1) * [(vx_b - 0), -(vy_b - 0);
                              (vy_b - 0),  (vx_b - 0)];
    masses.vp1_O = dRA_O * p1_A(1:2, :);
    
    % Acceleration Analysis for Link 1 Point Masses
    ddRA_O = (1 / robot.L1) * [(ax_b - 0), -(ay_b - 0);
                               (ay_b - 0),  (ax_b - 0)];
    masses.ap1_O = ddRA_O * p1_A(1:2, :);

    %=====================================================================%
    %                       LINK 2 POINT MASSES                           %
    %=====================================================================%
    
    % Non-rigid transformation to find points p2 in the CoM reference frame
    D2 = sqrt((1 / robot.m2) * robot.E2);
    p2 = D2 * q;
    
    % Center of mass coordinates of Link 2 relative to Joint B
    xcm2_B = robot.L2 / 2;
    ycm2_B = 0;
    
    % Homogeneous transformation from CoM to Joint B
    H2 = [1, 0, xcm2_B;
          0, 1, ycm2_B;
          0, 0, 1];
    p2_B = H2 * p2;
    
    % Transformation from Joint B to Inertial Frame (O)
    RB_O = (1 / robot.L2) * [(x_c - x_b), -(y_c - y_b);
                             (y_c - y_b),  (x_c - x_b)];
    HB_O = [RB_O, kin.B;
            0, 0, 1];
            
    masses.p2_O = HB_O * p2_B;
    
    % Velocity Analysis for Link 2 Point Masses
    dRB_O = (1 / robot.L2) * [(vx_c - vx_b), -(vy_c - vy_b);
                              (vy_c - vy_b),  (vx_c - vx_b)];
    masses.vp2_O = kin.vB + dRB_O * p2_B(1:2, :);
    
    % Acceleration Analysis for Link 2 Point Masses
    ddRB_O = (1 / robot.L2) * [(ax_c - ax_b), -(ay_c - ay_b);
                               (ay_c - ay_b),  (ax_c - ax_b)];
    masses.ap2_O = kin.aB + ddRB_O * p2_B(1:2, :);

end