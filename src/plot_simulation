function plot_simulation(t, results)
% PLOT_SIMULATION Visualizes the inverse dynamics simulation results.
%
% Generates standardized plots for the required joint torques computed via 
% Davies' method, isolating the visualization layer from the physics engine.
%
% Args:
%   t (1D array): Time vector of the simulation (s).
%   results (struct): Struct containing the computed dynamic outputs.
%       - t1 (1D array): Torque at joint A over time (N.m).
%       - t2 (1D array): Torque at joint B over time (N.m).

    % Plot Torque for Joint A
    figure('Name', 'Joint A Torque', 'NumberTitle', 'off');
    plot(t, results.t1, 'b', 'LineWidth', 3);
    title('Magnitude of Torque Applied in Coupling A (N.m) vs. Time (s)');
    xlabel('Time (s)');
    ylabel('Torque in Joint A (N.m)');
    grid on;
    % Tip: To compare with GIM, uncomment the lines below and load your CSV
    % hold on;
    % t1_GIM = table2array(readtable('t1.csv'))';
    % plot(t, t1_GIM, 'r*');
    % legend('Sim.', 'GIM');

    % Plot Torque for Joint B
    figure('Name', 'Joint B Torque', 'NumberTitle', 'off');
    plot(t, results.t2, 'b', 'LineWidth', 3);
    title('Magnitude of Torque Applied in Coupling B (N.m) vs. Time (s)');
    xlabel('Time (s)');
    ylabel('Torque in Joint B (N.m)');
    grid on;
    % Tip: To compare with GIM, uncomment the lines below and load your CSV
    % hold on;
    % t2_GIM = table2array(readtable('t2.csv'))';
    % plot(t, t2_GIM, 'r*');
    % legend('Sim.', 'GIM');
end