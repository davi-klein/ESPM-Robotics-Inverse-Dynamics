%=========================================================================%
% main_2R_simulation.m
% Refactored inverse dynamics simulation for the 2R planar robot.
%=========================================================================%
clc; close all; clear; format short eng;

% Adiciona a pasta de módulos ao path do MATLAB
addpath('src');

disp('Initializing 2R Robot Simulation...');

% 1. Initialization of Robot Parameters
robot = init_robot_params();

% 2. Simulation Parameters & Time Vector
n = 120;                  % Number of simulation points
s = 2.2;                  % Desired simulation time (s)
t = linspace(0, s, n+1);  % Time vector

% 3. Trajectory Generation (Circular Path & Inverse Kinematics)
th = 0 : pi/(n/2) : 2*pi;
rd = 0.5; x = 1; y = 0;
Px = x - rd * cos(th);
Py = y - rd * sin(th);

th2 = acos((Px.^2 + Py.^2 - robot.L1^2 - robot.L2^2) / (2 * robot.L1 * robot.L2));
th1 = atan2(Py, Px) - atan2((robot.L2 * sin(th2)), (robot.L1 + robot.L2 * cos(th2)));

% Input Joint Velocities and Accelerations (Defined in original research)
dth1 = linspace(pi/2, pi/2, n+1); 
ddth1 = zeros(1, n+1);            
dth2 = linspace(pi/2, pi/2, n+1); 
ddth2 = zeros(1, n+1);            

% 4. Data Structures to Store Results
results.t1 = zeros(1, n+1);
results.t2 = zeros(1, n+1);
results.pos = zeros(2, n+1);

%=========================================================================%
%                        MAIN DYNAMIC LOOP                                %
%=========================================================================%
disp('Computing Kinematics and Inverse Dynamics (Davies Method)...');

for i = 1:n+1
    % Step A: Kinematics
    kin = compute_kinematics(robot, th1(i), th2(i), dth1(i), dth2(i), ddth1(i), ddth2(i));
    
    % Step B: Equimomental Point Masses
    masses = compute_point_masses(robot, kin);
    
    % Step C: Inverse Dynamics 
    torques = compute_inverse_dynamics(robot, kin, masses);
    
    % Store Iteration Results
    results.pos(:, i) = kin.C;
    results.t1(i) = torques.joint_A;
    results.t2(i) = torques.joint_B;
end

%=========================================================================%
%                          VISUALIZATION                                  %
%=========================================================================%
disp('Generating plots...');
plot_simulation(t, results);

disp('Simulation finished successfully.');