function plot_simulation(t, results)
% PLOT_SIMULATION Visualizes the inverse dynamics simulation results.
%
% Dynamically generates standard plots for the required joint torques, 
% making it agnostic to the robot topology (2R, 3-RRR, 6-SPS).
%
% Args:
%   t (1D array): Time vector of the simulation (s).
%   results (struct): Struct containing dynamic outputs (e.g., t1, t2, t3).

    % Extract all fields (torques) from the results struct dynamically
    torque_fields = fieldnames(results);
    
    torque_fields = torque_fields(~strcmp(torque_fields, 'pos'));

    for idx = 1:length(torque_fields)
        field_name = torque_fields{idx};
        torque_data = results.(field_name);
        
        figure('Name', ['Joint ', num2str(idx), ' Torque'], 'NumberTitle', 'off');
        plot(t, torque_data, 'b', 'LineWidth', 3);
        title(['Magnitude of Torque Applied in Coupling ', num2str(idx), ' (N.m) vs. Time (s)']);
        xlabel('Time (s)');
        ylabel(['Torque in Joint ', num2str(idx), ' (N.m)']);
        grid on;
    end
end