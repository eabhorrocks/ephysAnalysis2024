function plot_wheel_with_states(wheelTime, wheelSpeed, stat_intervals, loco_intervals, ax)
% PLOT_WHEEL_WITH_STATES Plots wheel speed with colored background patches
%
% Usage:
%   plot_wheel_with_states(time, speed, stat_ints, loco_ints)
%   plot_wheel_with_states(time, speed, stat_ints, loco_ints, ax_handle)
%
% Inputs:
%   wheelTime:      Vector of time points
%   wheelSpeed:     Vector of speed values
%   stat_intervals: Nx2 matrix of [start, end] times for stationary
%   loco_intervals: Mx2 matrix of [start, end] times for locomotion
%   ax (optional):  Specific axes handle to plot into (default = gca)

    % 1. Handle optional axes argument
    if nargin < 5
        figure;
        ax = gca;
    end
    
    % Ensure we are adding to the specific axes
    axes(ax); 
    hold(ax, 'on');

    % 2. Determine Y-limits for the background patches
    % Add a 5% buffer so patches definitely cover the data
    y_range = max(wheelSpeed) - min(wheelSpeed);
    y_min = min(wheelSpeed) - (0.05 * y_range); 
    y_max = max(wheelSpeed) + (0.05 * y_range);
    
    % If signal is flat (0), create artificial height
    if y_min == y_max
        y_min = -1; y_max = 1;
    end

    % 3. Plot Stationary Intervals (Light Blue)
    if ~isempty(stat_intervals)
        % Create one patch object using NaN separators (faster than looping)
        % Strategy: Build X and Y vectors for all patches at once
        X_stat = [stat_intervals(:,1)'; stat_intervals(:,2)'; ...
                  stat_intervals(:,2)'; stat_intervals(:,1)'];
        Y_stat = repmat([y_min; y_min; y_max; y_max], 1, size(stat_intervals, 1));
        
        p_stat = patch(ax, X_stat, Y_stat, [0.8 0.8 1], ...
                      'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                      'DisplayName', 'Stationary');
    end

    % 4. Plot Locomotion Intervals (Light Red)
    if ~isempty(loco_intervals)
        X_loco = [loco_intervals(:,1)'; loco_intervals(:,2)'; ...
                  loco_intervals(:,2)'; loco_intervals(:,1)'];
        Y_loco = repmat([y_min; y_min; y_max; y_max], 1, size(loco_intervals, 1));
        
        p_loco = patch(ax, X_loco, Y_loco, [1 0.8 0.8], ...
                      'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                      'DisplayName', 'Locomotion');
    end

    % 5. Plot the Wheel Speed Trace (Black)
    plot(ax, wheelTime, wheelSpeed, 'k', 'LineWidth', 1.2, 'DisplayName', 'Speed');

    % 6. Formatting
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Speed (cm/s)');
    title(ax, 'Behavioral State Classification');
    axis(ax, 'tight');
    
    % Only show valid legend entries
    legend(ax, 'show', 'Location', 'best');
    
    hold(ax, 'off');
    
    xlim([min(wheelTime), max(wheelTime)])

end