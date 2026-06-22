%=========================================================================%
% main_3RRR_simulation.m
% Refactored inverse dynamics simulation for the 3-RRR parallel robot.
%=========================================================================%
clc; close all; clear; format short eng;

addpath('src');

disp('Initializing 3-RRR Robot Simulation...');

% 1. Initialization of Robot Parameters
robot = init_robot_3rrr_params();

% 2. Simulation Parameters & Time Vector
n = 120;                  
s = 2.2;                  
t = linspace(0, s, n+1);  

% 3. Trajectory Generation (Circular Path)
th = 0 : pi/(n/2) : 2*pi;
rd = 0.03; x = 0.25; y = 0.1083;
qx = x - rd * cos(th);
qy = y - rd * sin(th);

% Input Joint Velocities and Accelerations (Active Joints)
dth = repmat(pi/2, 3, n+1); % 3 joints, n+1 steps
ddth = zeros(3, n+1);

% Initial guess for the Newton-Raphson numerical solver
q0_guess = [0.125; 0.1; 0.375; 0.15; 0.25; 0.2];

% Data Structures to Store Results
results.t1 = zeros(1, n+1);
results.t2 = zeros(1, n+1);
results.t3 = zeros(1, n+1);

%=========================================================================%
%                        MAIN DYNAMIC LOOP                                %
%=========================================================================%
disp('Computing Kinematics and Inverse Dynamics (Davies Method)...');

for i = 1:n+1
    % Step A: Kinematics (Includes Newton-Raphson for closed loop)
    kin = compute_kinematics_3rrr(robot, qx(i), qy(i), dth(:,i), ddth(:,i), q0_guess);
    
    % Update the guess for the next iteration to speed up convergence
    q0_guess = kin.q0_next; 
    
    % Step B: Equimomental Point Masses (7 Bodies)
    masses = compute_point_masses_3rrr(robot, kin);
    
    % Step C: Inverse Dynamics (Davies Method for Parallel Robot)
    torques = compute_inverse_dynamics_3rrr(robot, kin, masses);
    
    % Store Iteration Results
    results.t1(i) = torques.joint_A1;
    results.t2(i) = torques.joint_A2;
    results.t3(i) = torques.joint_A3;
end

%=========================================================================%
%                          VISUALIZATION                                  %
%=========================================================================%
disp('Generating plots...');
plot_simulation(t, results);

disp('Simulation finished successfully.');
