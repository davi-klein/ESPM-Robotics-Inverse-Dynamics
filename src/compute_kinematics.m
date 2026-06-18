function kin = compute_kinematics(robot, th1, th2, dth1, dth2, ddth1, ddth2)
% COMPUTE_KINEMATICS Calculates forward and differential kinematics.
%
% Evaluates the Cartesian positions, velocities, and accelerations for 
% the 2R planar manipulator given the current joint states.
%
% Args:
%   robot (struct): The robot structural parameters.
%   th1, th2 (double): Joint angles in radians.
%   dth1, dth2 (double): Joint velocities in rad/s.
%   ddth1, ddth2 (double): Joint accelerations in rad/s^2.
%
% Returns:
%   kin (struct): Kinematic data including:
%       - B, C (2x1 double): Cartesian positions of joint B and end-effector C.
%       - vB, vC (2x1 double): Linear velocities of points B and C.
%       - aB, aC (2x1 double): Linear accelerations of points B and C.

    % Position Analysis
    kin.B = [robot.L1 * cos(th1); robot.L1 * sin(th1)];
    kin.C = [robot.L1 * cos(th1) + robot.L2 * cos(th1 + th2); 
             robot.L1 * sin(th1) + robot.L2 * sin(th1 + th2)];
    
    % Jacobian Matrix (Points-of-interest independent)
    Jj = [ robot.L1 * sin(th1), 0;
          -robot.L1 * cos(th1), 0; 
           robot.L1 * sin(th1) + robot.L2 * sin(th1 + th2), robot.L2 * sin(th1 + th2);
          -robot.L1 * cos(th1) - robot.L2 * cos(th1 + th2), -robot.L2 * cos(th1 + th2)]; 
          
    % Jacobian Matrix (Points-of-interest dependent)
    Jc = eye(4); 
    
    % Velocity Analysis
    v = Jc \ -(Jj * [dth1; dth2]);
    kin.vB = [v(1); v(2)];
    kin.vC = [v(3); v(4)];
    
    % Acceleration Analysis
    dJ = [0, 0, 0, 0, robot.L1 * cos(th1) * dth1, 0;
          0, 0, 0, 0, robot.L1 * sin(th1) * dth1, 0;
          0, 0, 0, 0, (robot.L1 * cos(th1) * dth1) + (robot.L2 * cos(th1 + th2) * (dth1 + dth2)), robot.L2 * cos(th1 + th2) * (dth1 + dth2);
          0, 0, 0, 0, (robot.L1 * sin(th1) * dth1) + (robot.L2 * sin(th1 + th2) * (dth1 + dth2)), robot.L2 * sin(th1 + th2) * (dth1 + dth2)];
          
    ac = Jc \ (-dJ * [v; dth1; dth2] - (Jj * [ddth1; ddth2]));
    kin.aB = [ac(1); ac(2)];
    kin.aC = [ac(3); ac(4)];
end