function robot = init_robot_6sps_params()
% INIT_ROBOT_6SPS_PARAMS Initializes properties of the 6-SPS parallel robot.
%
% Creates a structured parameter set for the spatial Stewart-Gough platform,
% organizing base geometry, moving platform dimensions, and the 3D inertial
% properties for the 13 rigid bodies (12 leg components + 1 platform).
%
% Returns:
%   robot (struct): Structural and inertial data including:
%       - g (double): Gravity.
%       - base (6x3 double): Coordinates of base joints A1 to A6.
%       - platform_pts (6x3 double): Coordinates of platform joints B1 to B6.
%       - legs (struct array): Mass, inertia, and pseudo-inertia for leg parts.
%       - platform (struct): Mass, inertia, and pseudo-inertia for the platform.

    % Environmental parameters
    robot.g = 9.807; % Note: Z-axis usually points up, gravity is positive down the math here

    % Fixed base joint coordinates (A1 to A6) [x, y, z]
    robot.base = [-2.120,  1.374, 0;
                  -2.380,  1.224, 0;
                  -2.380, -1.224, 0;
                  -2.120, -1.374, 0;
                   0.000, -0.150, 0;
                   0.000,  0.150, 0];

    % Moving platform local coordinates (B1_p to B6_p) [x, y, z]
    robot.platform_pts = [ 0.170,  0.595, -0.4;
                          -0.600,  0.150, -0.4;
                          -0.600, -0.150, -0.4;
                           0.170, -0.595, -0.4;
                           0.430, -0.445, -0.4;
                           0.430,  0.445, -0.4];

    % Leg Inertial Properties Initialization (Cylinders and Pistons)
    leg_mass = 0.1;
    leg_inertia = diag([0.00625, 0.00625, 0]); % Slender rod approximation
    
    for i = 1:6
        % Part 1 (e.g., Cylinder)
        robot.legs(i).part1.m = leg_mass;
        robot.legs(i).part1.I = leg_inertia;
        robot.legs(i).part1.xcm = 0; robot.legs(i).part1.ycm = 0; robot.legs(i).part1.zcm = 0;
        robot.legs(i).part1.E = compute_pseudo_inertia_3d(leg_inertia, leg_mass, 0, 0, 0);
        
        % Part 2 (e.g., Piston)
        robot.legs(i).part2.m = leg_mass;
        robot.legs(i).part2.I = leg_inertia;
        robot.legs(i).part2.xcm = 0; robot.legs(i).part2.ycm = 0; robot.legs(i).part2.zcm = 0;
        robot.legs(i).part2.E = compute_pseudo_inertia_3d(leg_inertia, leg_mass, 0, 0, 0);
    end

    % End-Effector (Moving Platform) Inertial Properties
    robot.platform.m = 1.5;
    robot.platform.I = diag([0.08, 0.08, 0.08]);
    robot.platform.xcm = 0;
    robot.platform.ycm = 0;
    robot.platform.zcm = 0;
    robot.platform.E = compute_pseudo_inertia_3d(robot.platform.I, robot.platform.m, 0, 0, 0);
end

function E = compute_pseudo_inertia_3d(I, m, xcm, ycm, zcm)
% COMPUTE_PSEUDO_INERTIA_3D Helper function for 3D equimomental system prep.
    E = [0.5*(-I(1,1)+I(2,2)+I(3,3)), -I(1,2), -I(1,3), m*xcm;
         -I(1,2), 0.5*(I(1,1)-I(2,2)+I(3,3)), -I(2,3), m*ycm;
         -I(1,3), -I(2,3), 0.5*(I(1,1)+I(2,2)-I(3,3)), m*zcm;
          m*xcm, m*ycm, m*zcm, m];
end