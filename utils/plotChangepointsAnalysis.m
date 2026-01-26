intervals1 = stationary_intervals;
intervals2 = locomotion_intervals;

% --- Define display properties
color1 = 'b'; % Blue for the first set
color2 = 'r'; % Red for the second set
transparency = 0.3; % Use a value between 0 (fully transparent) and 1 (fully opaque)
y_axis_height = [0 1]; % The vertical range the shaded regions will cover

% --- Create a new figure and hold it for multiple plots
figure;
hold on;


% --- Plot the first set of intervals
% We create a dummy handle 'p1' for the legend.
% It will be overwritten in the loop, and the last patch will be used for the legend entry.
p1 = patch(NaN, NaN, color1); % Dummy patch for the legend
for i = 1:size(intervals1, 1)
    start_time = intervals1(i, 1);
    end_time = intervals1(i, 2);
    
    % Define the vertices (corners) of the rectangle to be shaded
    % Format: [x1, x2, x3, x4], [y1, y2, y3, y4]
    vertices_x = [start_time, end_time,   end_time,   start_time];
    vertices_y = [y_axis_height(1), y_axis_height(1), y_axis_height(2), y_axis_height(2)];
    
    % Create the patch
    p1 = patch(vertices_x, vertices_y, color1, 'EdgeColor', 'none', 'FaceAlpha', transparency);
end

% --- Plot the second set of intervals
p2 = patch(NaN, NaN, color2); % Dummy patch for the legend
for i = 1:size(intervals2, 1)
    start_time = intervals2(i, 1);
    end_time = intervals2(i, 2);
    
    vertices_x = [start_time, end_time,   end_time,   start_time];
    vertices_y = [y_axis_height(1), y_axis_height(1), y_axis_height(2), y_axis_height(2)];
    
    p2 = patch(vertices_x, vertices_y, color2, 'EdgeColor', 'none', 'FaceAlpha', transparency);
end


% --- Finalize the plot
hold off; % Release the figure
box on; % Puts a box around the plot

% Set plot limits (optional, but good practice)
% This ensures all your intervals are visible.
all_intervals = [intervals1; intervals2];
xlim([min(all_intervals(:))-1, max(all_intervals(:))+1]);
ylim(y_axis_height);

% Add labels, title, and a legend
xlabel('Time or Value');
ylabel('State');
title('Shaded Regions for Two Sets of Intervals');
legend([p1, p2], {'Interval Set 1', 'Interval Set 2'}, 'Location', 'northwest');

% Remove the y-axis ticks for a cleaner look
set(gca, 'YTick', []);
hold on
plot(wheelTime, wheelSpeed, 'k')
