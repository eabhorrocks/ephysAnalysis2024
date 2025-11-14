function channel_row = depth_to_channel_row_np2(depth_um)
    % depth_to_channel_row Converts depth in µm to a channel row for a Neuropixels 2.0 probe.
    %
    %   This function calculates the corresponding channel row number from a
    %   physical depth along the probe shank. The calculation is based on the
    %   standard vertical pitch of Neuropixels 2.0 probes.
    %
    %   Syntax:
    %       channel_row = depth_to_channel_row(depth_um)
    %
    %   Input:
    %       depth_um (double): The depth along the probe in micrometers,
    %                          measured from the center of the lowest channel
    %                          row (row 0).
    %
    %   Output:
    %       channel_row (integer): The corresponding channel row number, rounded
    %                              down to the nearest integer.
    %
    %   Example:
    %       row = depth_to_channel_row(150);
    %       disp(row); % Displays 10

    % Neuropixels 2.0 probes have a vertical pitch of 15 µm between channel rows.
    VERTICAL_PITCH_UM = 15;

    % Input validation: ensure depth is not negative.
    if depth_um < 0
        error('Depth cannot be negative.');
    end

    % The channel row is determined by dividing the depth by the vertical
    % spacing. The floor function is used to round down to the nearest
    % complete row, which is equivalent to integer division for positive numbers.
    channel_row = floor(depth_um / VERTICAL_PITCH_UM);

end

% --- Example Usage ---
% To run this example, you can save the function above as 'depth_to_channel_row.m'
% and then run this script from the same directory or execute the lines in the
% MATLAB command window.

%{
    clc; % Clear command window

    % Example depths to test the conversion
    example_depths = [0, 14, 15, 29, 30, 150, 765.5];

    fprintf('--- Neuropixels 2.0 Depth to Channel Row Conversion ---\n');
    fprintf('Assuming a vertical pitch of 15 µm between rows.\n\n');

    for i = 1:length(example_depths)
        depth = example_depths(i);
        try
            row = depth_to_channel_row(depth);
            fprintf('Depth: %7.1f µm  ->  Channel Row: %d\n', depth, row);
        catch ME
            fprintf('Error processing depth %.1f µm: %s\n', depth, ME.message);
        end
    end

    % Example of handling invalid input
    fprintf('\nTesting with negative depth:\n');
    try
        depth_to_channel_row(-10);
    catch ME
        fprintf('Caught expected error: %s\n', ME.message);
    end
%}
