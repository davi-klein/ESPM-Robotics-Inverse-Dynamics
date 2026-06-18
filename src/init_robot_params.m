function robot = init_robot_params()
% INIT_ROBOT_PARAMS Initializes the physical properties of the 2R robot.
%
% This function sets up the structural parameters, masses, and inertia
% tensors required for the inverse dynamics simulation using the Davies
% method and equimomental systems.
%
% Returns:
%   robot (struct): A structured variable containing:
%       - g (double): Gravitational acceleration (m/s^2).
%       - L1, L2 (double): Link lengths (m).
%       - m1, m2 (double): Link masses (kg).
%       - I1, I2 (double): 3x3 inertia tensors.
%       - E1, E2 (double): 3x3 pseudo-inertia matrices.

    % Environmental parameters
    robot.g = -9.81;
    
    % Link 1 properties
    robot.L1 = 1.0;
    robot.m1 = 1.0;
    robot.I1 = diag([0.08333333, 0.08333333, 0.08333333]);
    robot.xcm1 = 0;
    robot.ycm1 = 0;
    
    % Link 2 properties
    robot.L2 = 1.0;
    robot.m2 = 1.0;
    robot.I2 = diag([0.08333333, 0.08333333, 0.08333333]);
    robot.xcm2 = 0;
    robot.ycm2 = 0;
    
    % Pseudo-inertia matrices (Equimomental system preparation)
    robot.E1 = compute_pseudo_inertia(robot.I1, robot.m1, robot.xcm1, robot.ycm1);
    robot.E2 = compute_pseudo_inertia(robot.I2, robot.m2, robot.xcm2, robot.ycm2);
end

function E = compute_pseudo_inertia(I, m, xcm, ycm)
% COMPUTE_PSEUDO_INERTIA Helper function to calculate pseudo-inertia matrix.
    E = [0.5*(-I(1,1)+I(2,2)+I(3,3)), -I(1,2), m*xcm;
         -I(1,2), 0.5*(I(1,1)-I(2,2)+I(3,3)),  m*ycm;
          m*xcm, m*ycm, m];
end