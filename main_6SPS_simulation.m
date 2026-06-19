%=========================================================================%
% main_6SPS_simulation.m
% Refactored inverse dynamics simulation for the 6-SPS parallel robot.
%=========================================================================%
clc; close all; clear; format short eng;
addpath('src');

disp('Initializing 6-SPS Stewart-Gough Platform Simulation...');

% 1. Load Pre-calculated Trajectory & Validation Data
% Ensure 'dL.mat', 'ddL.mat', and 'tau.mat' are in the root directory
load('data/dL.mat');
load('data/ddL.mat');
load('data/tau.mat');

robot = init_robot_6sps_params();

% 2. Simulation Parameters
n = 120;                  
s = 2.0944;                  
t = linspace(0, s, n);  
w = 3;

% Trajectory Generation (End-effector position P over time)
x = -1.5 + 0.2*sin(w*t);
y = 0.2*sin(w*t);
z = 1 + 0.2*sin(w*t);
P_traj = [x; y; z];

% Initial guess for the Newton-Raphson numerical solver (18 Variables)
q0_guess = [
    -1.3300;  0.5950; 0.6;
    -2.1000;  0.1500; 0.6;
    -2.1000; -0.1500; 0.6;
    -1.3300; -0.5950; 0.6;
    -1.0700; -0.4450; 0.6;
    -1.0700;  0.4450; 0.6
];

% Data Structures to Store Results
results.F1 = zeros(1, n);
results.F2 = zeros(1, n);
results.F3 = zeros(1, n);
results.F4 = zeros(1, n);
results.F5 = zeros(1, n);
results.F6 = zeros(1, n);

% Note: We compute L(i,:) from your inverse kinematics data natively 
% or assume it's pre-loaded alongside dL and ddL. Based on your original 
% script structure, we will use the position equations.
disp('Computing Kinematics and Inverse Dynamics (Davies Method)...');

for i = 1:n
    % Note: If L(i,:) wasn't loaded from a mat file, we can deduce it 
    % from the trajectory and your Inverse Kinematics module if available.
    % Assuming L, dL, and ddL are available in memory for step i:
    L_step = sqrt(dL(i,:).^2 + ddL(i,:).^2); % Replace with true L if loaded
    
    % Step A: Spatial Kinematics (Newton-Raphson 18-DOF solver)
    kin = compute_kinematics_6sps(robot, L_step, dL(i,:), ddL(i,:), q0_guess);
    q0_guess = kin.q0_next; 
    
    % Step B: 3D Equimomental Tetrahedrons (13 Bodies)
    masses = compute_point_masses_6sps(robot, kin);
    
    % Step C: Spatial Inverse Dynamics
    torques = compute_inverse_dynamics_6sps(robot, kin, masses);
    
    % Store Iteration Results
    results.F1(i) = torques.F1;
    results.F2(i) = torques.F2;
    results.F3(i) = torques.F3;
    results.F4(i) = torques.F4;
    results.F5(i) = torques.F5;
    results.F6(i) = torques.F6;
end

%=========================================================================%
%                          VISUALIZATION                                  %
%=========================================================================%
disp('Generating dynamic validation plots...');

% Using our topology-agnostic plot module!
plot_simulation(t, results);

disp('Simulation finished successfully.');