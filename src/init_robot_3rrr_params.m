function robot = init_robot_3rrr_params()
% INIT_ROBOT_3RRR_PARAMS Initializes properties of the 3-RRR parallel robot.
%
% This function creates a structured parameter set for the 3-RRR manipulator,
% organizing base geometry, moving platform dimensions, and the inertial
% properties for the 6 leg links and the end-effector.
%
% Returns:
%   robot (struct): Structural and inertial data including:
%       - g (double): Gravity.
%       - base (3x2 double): Coordinates of base joints A1, A2, A3.
%       - d (1x3 double): Moving platform joint distances d1, d2, d3.
%       - L (1x6 double): Lengths of the 6 leg links.
%       - legs (struct array): Mass, inertia, and pseudo-inertia for legs.
%       - platform (struct): Mass, inertia, and pseudo-inertia for platform.

    % Environmental parameters
    robot.g = -9.81;

    % Fixed base joint coordinates (A1, A2, A3)
    robot.base = [0, 0; 
                  0.5, 0; 
                  0.25, 0.2886751];
                  
    % Moving platform geometry
    robot.d = [1/8, 1/8, 1/8];
    
    % Link lengths (L1 through L6)
    robot.L = (1/6) * ones(1, 6);

    % Leg Inertial Properties Initialization
    base_inertia = diag([0.002314734, 0.002314734, 0.002314734]);
    
    for i = 1:6
        robot.legs(i).m = 1.0;
        robot.legs(i).I = base_inertia;
        robot.legs(i).xcm = 0; % CoM relative to local frame
        robot.legs(i).ycm = 0;
        robot.legs(i).xcm_joint = robot.L(i) / 2; % CoM relative to previous joint
        robot.legs(i).ycm_joint = 0;
        
        % Pre-calculate Pseudo-Inertia
        robot.legs(i).E = compute_pseudo_inertia(robot.legs(i).I, robot.legs(i).m, ...
                                                 robot.legs(i).xcm, robot.legs(i).ycm);
    end

    % End-Effector (Moving Platform) Inertial Properties
    robot.platform.m = 1.0;
    robot.platform.I = diag([0.001302093, 0.001302093, 0.001302093]);
    robot.platform.xcm = 0;
    robot.platform.ycm = 0;
    robot.platform.E = compute_pseudo_inertia(robot.platform.I, robot.platform.m, ...
                                              robot.platform.xcm, robot.platform.ycm);
end

function E = compute_pseudo_inertia(I, m, xcm, ycm)
% COMPUTE_PSEUDO_INERTIA Helper function for equimomental system prep.
    E = [0.5*(-I(1,1)+I(2,2)+I(3,3)), -I(1,2), m*xcm;
         -I(1,2), 0.5*(I(1,1)-I(2,2)+I(3,3)),  m*ycm;
          m*xcm, m*ycm, m];
end