%% Review batch LFP Sparse Noise

dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';
dataFiles = dir(fullfile(dataDir, 'M26001_20260122*_LFP_RFs_wLatency*')); % 
tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'lfp');
    % load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name),  'trials'); % just trials
    session(ifile).lfp = lfp;
    % session(ifile).trials = trials;
end
toc

lfp = cat(1,session.lfp);


%% Plot Azimuth as a Function of Channel Depth 

% Parameters
p_val_thresh = 0.05;
nonSigAlpha = 0.3;               % 0 = Invisible (White), 0.2 = Faint, 1 = Fully Visible
x_lims = [-120, 120];         % Visual Azimuth limits
unique_shanks = [0, 1, 2, 3]; % NP 2.0 usually has 4 shanks
nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 600]);
t = tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through each Shank
for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter channels belonging to this shank
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    % Check if we have data for this shank
    if isempty(shank_chans)
        nexttile; title(sprintf('Shank %d (Empty)', target_shank)); axis off;
        continue;
    end
    
    % 2. Sort these channels by Depth (Tip to Top)
    depths = [shank_chans.y];
    [sorted_depths, sort_idx] = sort(depths, 'ascend');
    
    % 3. Build the Image Stack & Significance Mask
    sample_rf = shank_chans(1).rf_map_contrast_scaleAndSmoothed;
    [~, map_width] = size(sample_rf);
    
    nChans = length(shank_chans);
    azimuth_stack = zeros(nChans, map_width);
    sig_mask_stack = zeros(nChans, 1); % Column vector for row transparency
    
    for k = 1:nChans
        % Get index in the shank_chans array
        orig_idx = sort_idx(k);
        
        % Get Map
        map = shank_chans(orig_idx).rf_map_contrast_scaleAndSmoothed;
        
        % Collapse Elevation (Min/Sink)
        azimuth_profile = min(map, [], 1);
        azimuth_stack(k, :) = azimuth_profile;
        
        % Check Significance Field (Safety check in case field missing)
        if isfield(shank_chans(orig_idx), 'global_p_value')
            sig_mask_stack(k) = shank_chans(orig_idx).global_p_value<=p_val_thresh;
        else
            sig_mask_stack(k) = 1; % Default to visible if no stats run
        end
    end
    
    % 4. Visualization
    nexttile;
    
    y_lims = [min(sorted_depths), max(sorted_depths)];
    
    % Prepare Plot Data (Invert so sink is positive)
    plot_data = azimuth_stack * -1;
    % Clamp negative values (sources) to 0 so colormap starts at 0
    plot_data(plot_data < 0) = 0;
    
    % Create Image Object
    im = imagesc(x_lims, y_lims, plot_data);
    
   % --- ALPHA MAPPING ---
    % 1. Create a vector of alpha values based on significance
    alpha_vals = sig_mask_stack; % Currently 0s and 1s
    
    % 2. Replace 0s (Non-Sig) with your parameter value
    alpha_vals(alpha_vals == 0) = nonSigAlpha;
    
    % 3. Replicate across width to make the full Alpha Map
    alpha_map = repmat(alpha_vals, 1, map_width);
    
    % 4. Apply
    set(im, 'AlphaData', alpha_map);
    % ---------------------
    
    set(gca, 'YDir', 'normal'); 
    colormap(parula);
    title(sprintf('Shank %d', target_shank));
    xlabel('Azimuth (deg)');
    
    if i == 1
        ylabel('Distance from tip (\mu m)');
    else
        % Remove Y-labels for inner plots to save space
        set(gca, 'YTickLabel', []);
    end
    
    % 5. Trend Line (Only for significant channels)
    hold on;
    [~, max_col_idx] = max(plot_data, [], 2);
    
    deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
    peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_px;
    
    % Filter: Only plot points where the channel is significant
    sig_indices = find(sig_mask_stack);
    
    % Plot peaks as dots (better than line when data is sparse)
    if ~isempty(sig_indices)
        scatter(peak_azimuths(sig_indices), sorted_depths(sig_indices), ...
                10, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
    hold off;
    xline(-120:20:120, 'w-', 'LineWidth', 1, 'Alpha', 0.5);

    
end

% Global Elements
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Sink Strength (Norm)';
sgtitle('Retinotopy by Shank (Significant Channels Only)');

%% Plot Peak Latency vs Depth (Significant Channels Only)

% --- Parameters ---
unique_shanks = [0, 1, 2, 3];
latency_lims  = [0, 150];     % X-axis limits in milliseconds (e.g., 0 to 150ms)
dot_size      = 25;           % Size of scatter points
% ------------------

nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 600]);
t = tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter Channels: Belong to Shank AND Significant
    %    We assume lfp is the struct array for the current probe
    %    (If lfp is 6x384, select one row first, e.g., lfp_row = lfp(1,:))
    
    mask = [lfp.shanks] == target_shank & [lfp.is_visual_channel] == 1;
    sig_chans = lfp(mask);
    
    nexttile;
    
    if isempty(sig_chans)
        title(sprintf('Shank %d (No Sig Data)', target_shank)); 
        axis off;
        continue;
    end
    
    % 2. Extract Data
    y_depths = [sig_chans.y];
    latencies_s = [sig_chans.peak_latency];
    
    % Convert to milliseconds for nicer plotting
    latencies_ms = latencies_s * 1000;
    
    % 3. Scatter Plot
    % X = Latency, Y = Depth
    scatter(latencies_ms, y_depths, dot_size, 'filled', ...
            'MarkerFaceColor', [0.2 0.2 0.2], ... % Dark Grey dots
            'MarkerFaceAlpha', 0.6);
        
    % 4. Formatting
    grid on;
    xlabel('Latency (ms)');
    xlim(latency_lims);
    ylim([min([lfp.y]), max([lfp.y])]); % Consistent depth axis across plots
    
    title(sprintf('Shank %d', target_shank));
    set(gca, 'YDir', 'normal'); % Deepest (0) at bottom
    
    % Add "Stimulus Onset" line
    xline(0, 'r--', 'Stimulus', 'LabelVerticalAlignment', 'bottom');
    
    % Only label Y-axis on the first plot
    if i == 1
        ylabel('Depth (\mu m)');
    else
        set(gca, 'YTickLabel', []);
    end
    
    % 5. Optional: Fit a simple trend line (Polynomial) to see the profile
    %    (Laminar profiles are often U-shaped, so a quadratic fit might work)
    % if length(latencies_ms) > 5
    %     [p, S] = polyfit(y_depths, latencies_ms, 2); 
    %     y_fit = linspace(min(y_depths), max(y_depths), 100);
    %     x_fit = polyval(p, y_fit);
    % 
    %     hold on;
    %     plot(x_fit, y_fit, 'b-', 'LineWidth', 1.5); % Blue trend line
    %     hold off;
    % end
end

sgtitle('Response Latency vs. Depth (Significant Channels)');