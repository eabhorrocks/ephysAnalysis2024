%% Review batch LFP Sparse Noise

dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';
dataFiles = dir(fullfile(dataDir, 'M26018_20260621*_LFP_RFs_wLatency*')); % 


% fields_to_drop = {'peak_latency', 'rf_latency_map'};

tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'lfp_save');
    % load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name),  'trials'); % just trials
    session(ifile).lfp = lfp_save;
    % session(ifile).trials = trials;

    % current_fields = fieldnames(session(ifile).lfp);
    % 
    % % Only attempt to remove the fields that are actually present
    % present_to_drop = intersect(fields_to_drop, current_fields);
    % 
    % if ~isempty(present_to_drop)
    %     session(ifile).lfp = rmfield(session(ifile).lfp, present_to_drop);
    % end
end
toc

lfp = cat(2,session.lfp);


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
    xlim([-125 125])
    ax=gca;
    ax.XTick = -120:20:120;
    
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


%% Complete CCF Atlas Overlay & Trajectory Reconstruction (Corrected Geometry)

% =========================================================================
% 1. CONFIGURATION & INSERTION PARAMETERS
% =========================================================================
% Target Stereotaxic Coordinates (Bregma = 0)
insertion_AP = -3400;   
insertion_ML = 949;     % Right Hemisphere (+ is right)
insertion_DV = 124;     
insertion_depth = 4900; % Total distance driven into brain (um)

% CCF Atlas Bregma Offsets (Absolute 10um Pixel Space)
ccf_bregma_AP = 5400; 
ccf_bregma_ML = 5700;   % Midline is at 5700 um
ccf_bregma_DV = 0;      % Top of the CCF volume

% Plotting Flags
plot_layers = false;    % FALSE = Fused Parent Regions. TRUE = Individual Layers.

% =========================================================================
% 2. COORDINATE TRANSFORMATION (PROBE -> CCF SPACE)
% =========================================================================
fprintf('Transforming probe coordinates to CCF space...\n');
y_coords = [lfp.y];
x_coords = [lfp.x];

% Anchor the First Shank to Origin
x_coords = x_coords - min(x_coords);

% Distance down the probe shaft from the surface
d_surface = insertion_depth - y_coords; 

% 45-Degree Lateral Tilt Geometry
tilt_deg = 45; 

% 1. Shaft Trajectory: Probe points Down (+DV) and Right/Lateral (+ML)
delta_ML = d_surface .* sind(tilt_deg);
delta_DV = d_surface .* cosd(tilt_deg);

% 2. Shank Spacing: Spread the 4 shanks perpendicularly across the coronal plane
% To separate them visually, they stagger Up (-DV) and Right (+ML)
offset_ML = x_coords .* cosd(tilt_deg);
offset_DV = -x_coords .* sind(tilt_deg);

% Combine to get True Stereotaxic Coordinates
stereo_ML = insertion_ML + delta_ML + offset_ML;
stereo_DV = insertion_DV + delta_DV + offset_DV;

% Shift to absolute CCF Atlas Space
ccf_ML = ccf_bregma_ML + stereo_ML;
ccf_DV = ccf_bregma_DV + stereo_DV;

for i = 1:numel(lfp)
    lfp(i).ML = ccf_ML(i);
    lfp(i).DV = ccf_DV(i);
end

% =========================================================================
% 3. LOAD CCF VOLUMES & EXTRACT CORONAL SLICE
% =========================================================================
fprintf('Checking CCF Volumes...\n');

% Calculate the CCF Slice Index based on your AP coordinate
ap_slice_idx = round((-insertion_AP + ccf_bregma_AP) / 10); 

% Only load the massive CCF files if they are NOT already in the workspace
if ~exist('tv', 'var') || ~exist('anno', 'var') || ~exist('st', 'var')
    fprintf('  -> Loading Allen CCF Volumes from disk (this takes a moment)...\n');
    tv = readNPY('template_volume_10um.npy'); 
    anno = readNPY('annotation_volume_10um_by_index.npy');
    st = readtable('structure_tree_safe_2017.csv');
else
    fprintf('  -> CCF Volumes found in workspace. Skipping load step!\n');
end

% Always extract the 2D Coronal Slices (in case the AP coordinate was changed)
coronal_slice = squeeze(tv(ap_slice_idx, :, :));
anno_slice = squeeze(anno(ap_slice_idx, :, :));

% Isolate the Right Hemisphere 
right_hemi_mask = false(size(anno_slice));
right_hemi_mask(:, 570:end) = true; % Pixel 570 is the CCF midline
anno_slice_right = anno_slice;
anno_slice_right(~right_hemi_mask) = 0;

% =========================================================================
% 4. EXTRACT BOUNDARIES FOR VISUAL & RETROSPLENIAL AREAS
% =========================================================================
fprintf('Extracting anatomical boundaries...\n');
unique_indices = unique(anno_slice_right);
unique_indices(unique_indices == 0) = []; 

target_areas = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(unique_indices)
    current_idx = unique_indices(i);
    
    if current_idx > height(st) || current_idx < 1
        continue; 
    end
    acronym = st.acronym{current_idx};
    
    % Grab all Visual Cortex (exclude Visceral) AND all Retrosplenial Cortex
    isVisual = startsWith(acronym, 'VIS') && ~startsWith(acronym, 'VISC');
    isRSP    = startsWith(acronym, 'RSP');
    
    if isVisual || isRSP
        if plot_layers
            plot_acronym = acronym; 
        else
            plot_acronym = regexprep(acronym, '[0-9/].*', '');
        end
        
        layer_mask = (anno_slice_right == current_idx);
        
        if isKey(target_areas, plot_acronym)
            target_areas(plot_acronym) = target_areas(plot_acronym) | layer_mask;
        else
            target_areas(plot_acronym) = layer_mask;
        end
    end
end

% =========================================================================
% 5. VISUALIZATION
% =========================================================================
fprintf('Generating plot...\n');
figure('Color', 'w', 'Position', [200, 100, 1000, 800]);

% A. Draw faint grey morphological background
imagesc(coronal_slice); colormap(gray); hold on;

% B. Draw faint white boundaries for ALL regions
all_boundaries = edge(double(anno_slice_right), 'Sobel');
[all_y, all_x] = find(all_boundaries);
scatter(all_x, all_y, 1, [1 1 1], 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1); 

% C. Draw Bold Boundaries & Acronyms
area_keys = keys(target_areas);
colors = jet(length(area_keys)); 

for k = 1:length(area_keys)
    area_name = area_keys{k};
    area_mask = target_areas(area_name);
    
    boundaries = edge(double(area_mask), 'Sobel');
    [edge_y, edge_x] = find(boundaries);
    scatter(edge_x, edge_y, 6, colors(k,:), 'filled'); 
    
    [y_coords, x_coords] = find(area_mask);
    if ~isempty(x_coords)
        centroid_x = mean(x_coords);
        centroid_y = mean(y_coords);
        
        text(centroid_x, centroid_y, area_name, ...
            'Color', colors(k,:), 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Interpreter', 'none');
    end
end

% D. Plot the 4-Shank Probe Trajectory
scatter([lfp.ML] / 10, [lfp.DV] / 10, 3, [lfp.shanks], 'filled', 'r');

% E. Formatting
axis equal; 
set(gca, 'YDir', 'reverse'); 
xlim([400, 1140]); 
ylim([0, 800]);

if plot_layers
    title_str = sprintf('Visual & Retrosplenial Layers (AP \\approx %d \\mu m)', insertion_AP);
else
    title_str = sprintf('Visual & Retrosplenial Areas (AP \\approx %d \\mu m)', insertion_AP);
end
title(title_str, 'FontSize', 14);
xlabel('Medial-Lateral (CCF Pixels)'); 
ylabel('Dorsal-Ventral (CCF Pixels)');

fprintf('Done.\n');

% Review Batch LFP Sparse Noise with Anatomical Axis (Option 1)

% =========================================================================
% 1. ANATOMICAL LOOKUP (Map channels to CCF Brain Areas)
% =========================================================================
fprintf('Extracting anatomical areas for each channel...\n');

% Ensure CCF Volumes are loaded (Requires them to be in the workspace)
if ~exist('anno', 'var') || ~exist('st', 'var')
    error('CCF Volumes not found. Please run the CCF loading script first!');
end

% Extract the 2D Coronal Slice of the Annotation Volume
ap_slice_idx = round((-insertion_AP + ccf_bregma_AP) / 10); 
anno_slice = squeeze(anno(ap_slice_idx, :, :));

% Assign brain area acronyms to the lfp struct
for i = 1:numel(lfp)
    % Convert coordinates (um) to CCF pixel indices (10um resolution)
    ml_idx = round(lfp(i).ML / 10);
    dv_idx = round(lfp(i).DV / 10);
    
    % Bounds check to prevent crashing if a channel is outside the volume
    if dv_idx > 0 && dv_idx <= size(anno_slice, 1) && ...
       ml_idx > 0 && ml_idx <= size(anno_slice, 2)
   
        % Look up the ID in the CCF volume
        pixel_id = anno_slice(dv_idx, ml_idx);
        
        if pixel_id > 0 && pixel_id <= height(st)
            acr = st.acronym{pixel_id};
            
            % If it's Visual or RSP cortex, strip the layer numbers (e.g. VISp1 -> VISp)
            if startsWith(acr, 'VIS') || startsWith(acr, 'RSP')
                acr = regexprep(acr, '[0-9/].*', '');
            end
            lfp(i).area = acr;
        else
            lfp(i).area = 'void';
        end
    else
        lfp(i).area = 'void';
    end
end

% Create a color map for the unique areas found across the probe
unique_areas = unique({lfp.area});
area_colors = containers.Map(unique_areas, mat2cell(lines(length(unique_areas)), ones(length(unique_areas),1), 3));

% =========================================================================
% 2. PLOT RETINOTOPY WITH ANATOMICAL AXIS
% =========================================================================
fprintf('Generating Functional-Anatomical Plot...\n');

% Parameters
p_val_thresh = 0.05;
nonSigAlpha = 0.3;               
x_lims = [-120, 120];         
unique_shanks = [0, 1, 2, 3]; 
nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 700]);
tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter channels belonging to this shank
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    if isempty(shank_chans)
        nexttile; title(sprintf('Shank %d (Empty)', target_shank)); axis off;
        continue;
    end
    
    % 2. Sort these channels by Depth (Tip to Top)
    depths = [shank_chans.y];
    [sorted_depths, sort_idx] = sort(depths, 'ascend');
    sorted_chans = shank_chans(sort_idx);
    
    % 3. Build the Image Stack & Significance Mask
    sample_rf = sorted_chans(1).rf_map_contrast_scaleAndSmoothed;
    [~, map_width] = size(sample_rf);
    
    nChans = length(sorted_chans);
    azimuth_stack = zeros(nChans, map_width);
    sig_mask_stack = zeros(nChans, 1); 
    
    for k = 1:nChans
        map = sorted_chans(k).rf_map_contrast_scaleAndSmoothed;
        azimuth_stack(k, :) = min(map, [], 1); % Collapse Elevation (Sink)
        
        if isfield(sorted_chans(k), 'global_p_value')
            sig_mask_stack(k) = sorted_chans(k).global_p_value <= p_val_thresh;
        else
            sig_mask_stack(k) = 1; 
        end
    end
    
    % Prepare Plot Data (Invert so sink is positive)
    plot_data = azimuth_stack * -1;
    plot_data(plot_data < 0) = 0; % Clamp negative values
    
    % 4. Visualization
    nexttile;
    y_lims = [min(sorted_depths), max(sorted_depths)];
    
    % Draw the functional heatmap
    im = imagesc(x_lims, y_lims, plot_data);
    
    % Apply Alpha Mask
    alpha_vals = sig_mask_stack;
    alpha_vals(alpha_vals == 0) = nonSigAlpha;
    alpha_map = repmat(alpha_vals, 1, map_width);
    set(im, 'AlphaData', alpha_map);
    
    set(gca, 'YDir', 'normal'); 
    colormap(gca, parula);
    hold on;
    
    % =====================================================================
    % 5. ADD ANATOMICAL AXIS (The colored bar on the left)
    % =====================================================================
    % Define the physical location of the colored bar
    bar_width = 15;
    bar_left = x_lims(1) - bar_width - 5; 
    bar_right = x_lims(1) - 5;
    
    % Extract areas for this specific shank
    shank_areas = {sorted_chans.area};
    
    % Find boundaries where the area changes
    boundary_idx = [1, find(~strcmp(shank_areas(1:end-1), shank_areas(2:end))) + 1, length(shank_areas)+1];
    
    for b = 1:(length(boundary_idx)-1)
        start_idx = boundary_idx(b);
        end_idx = boundary_idx(b+1) - 1;
        
        area_name = shank_areas{start_idx};
        area_color = area_colors(area_name);
        
        % Calculate Y limits for this anatomical block
        y_start = sorted_depths(start_idx);
        y_end = sorted_depths(end_idx);
        
        % Ignore 'void' areas or very tiny slivers to keep the plot clean
        if ~strcmp(area_name, 'void')
            % Draw the colored anatomical block
            fill([bar_left, bar_right, bar_right, bar_left], ...
                 [y_start, y_start, y_end, y_end], ...
                 area_color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
             
            % Draw horizontal dashed line across the heatmap
            yline(y_start, '--', 'Color', [1 1 1 0.5], 'LineWidth', 1);
            
            % Add Text Label inside the anatomical block
            mid_y = (y_start + y_end) / 2;
            text(bar_left - 5, mid_y, area_name, 'Color', area_color, ...
                 'FontSize', 9, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end
    end
    
    % =====================================================================
    % 6. Plot RF Peaks
    [~, max_col_idx] = max(plot_data, [], 2);
    deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
    peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_px;
    
    sig_indices = find(sig_mask_stack);
    if ~isempty(sig_indices)
        scatter(peak_azimuths(sig_indices), sorted_depths(sig_indices), ...
                10, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
    
    % Gridlines & Formatting
    xline(-120:20:120, 'w-', 'LineWidth', 1, 'Alpha', 0.3);

    % HORIZONTAL LINES AT QUARTILE DEPTHS ---
    y_span = y_lims(2) - y_lims(1);
    y_quartiles = y_lims(1) + y_span .* [0.25, 0.5, 0.75];
    yline(y_quartiles, 'r:', 'LineWidth', 1.5);
    
    % Expand X-limits to fit the new anatomical bar and text!
    xlim([bar_left - 45, x_lims(2)]); 
    ax = gca;
    ax.XTick = -120:40:120;
    
    title(sprintf('Shank %d', target_shank));
    xlabel('Azimuth (deg)');
    
    if i == 1
        ylabel('Distance from tip (\mu m)');
    else
        set(gca, 'YTickLabel', []);
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Sink Strength (Norm)';
sgtitle('Retinotopy by Shank with Anatomical Annotations', 'FontSize', 16, 'FontWeight', 'bold');

%% pinpoint version

%% Review Batch LFP Sparse Noise with Low-Level Pinpoint Parser

% =========================================================================
% 1. ANATOMICAL LOOKUP (Hardware ID to Spatial Depth Decoding)
% =========================================================================
fprintf('Parsing serialized Pinpoint text string...\n');

plot_layers = false; 

% Open file dialog to select the Pinpoint text file
[filename, filepath] = uigetfile('*.txt', 'Select your Pinpoint Channels Text File');
if isequal(filename,0)
    error('User cancelled file selection.');
else
    raw_text = fileread(fullfile(filepath, filename));
end

% Clean up the JSON tracking wrappers if present
raw_text = regexprep(raw_text, '^\["', '');
raw_text = regexprep(raw_text, '"\]$', '');

% Drop the session hash (e.g., "a26abd7a:")
string_parts = split(raw_text, ':');
if length(string_parts) > 1
    data_stream = string_parts{2};
else
    data_stream = string_parts{1};
end

% Split the stream into individual channel entries
channel_entries = split(data_stream, ';');

% --- RECONSTRUCT PINPOINT'S PHYSICAL MAP ---
pinpoint_map = struct('shank', [], 'y', [], 'area', {});
valid_count = 0;

for ch = 1:length(channel_entries)
    current_entry = channel_entries{ch};
    if isempty(strtrim(current_entry)), continue; end
    
    elements = split(current_entry, ',');
    
    if length(elements) >= 3
        valid_count = valid_count + 1;
        absolute_id = str2double(elements{1});
        acr = elements{3};
        
        % Clean up layer numbers if requested
        if ~plot_layers && (startsWith(acr, 'VIS') || startsWith(acr, 'RSP'))
            acr = regexprep(acr, '[0-9/].*', ''); 
        end
        
        % --- NEUROPIXELS 2.0 GEOMETRY MATH ---
        % Shank 0: 0-1279 | Shank 1: 1280-2559 | Shank 2: 2560-3839 | Shank 3: 3840-5119
        pp_shank = floor(absolute_id / 1280);
        
        % Local ID on that specific shank (0 to 1279)
        local_id = mod(absolute_id, 1280);
        
        % Vertical Depth: 2 contacts per row, rows are 15 um apart
        pp_y = floor(local_id / 2) * 15;
        
        pinpoint_map(valid_count).shank = pp_shank;
        pinpoint_map(valid_count).y = pp_y;
        pinpoint_map(valid_count).area = acr;
    end
end
fprintf('  -> Reconstructed spatial map for %d Pinpoint contacts.\n', valid_count);

% --- MAP TO LFP SUBSET (Spatial Nearest-Neighbor) ---
for i = 1:numel(lfp)
    chan_shank = lfp(i).shanks;
    chan_y     = lfp(i).y;
    
    % Filter the Pinpoint reference map to only look at the current shank
    shank_mask = [pinpoint_map.shank] == chan_shank;
    ref_shank_data = pinpoint_map(shank_mask);
    
    if isempty(ref_shank_data)
        lfp(i).area = 'void';
        continue;
    end
    
    % Find the Pinpoint reference contact physically closest to this LFP channel's depth
    [~, min_idx] = min(abs([ref_shank_data.y] - chan_y));
    
    lfp(i).area = ref_shank_data(min_idx).area;
end

% Generate dynamic plot colors for the unique areas detected
unique_areas = unique({lfp.area});
area_colors = containers.Map(unique_areas, mat2cell(lines(length(unique_areas)), ones(length(unique_areas),1), 3));

fprintf('Anatomical mapping complete.\n');


% =========================================================================
% 2. PLOT RETINOTOPY WITH ANATOMICAL AXIS
% =========================================================================
fprintf('Generating Functional-Anatomical Plot...\n');
% Parameters
p_val_thresh = 0.05;
nonSigAlpha = 0.3;               
x_lims = [-120, 120];         
unique_shanks = [0, 1, 2, 3]; 
nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 700]);
tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter channels belonging to this shank
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    if isempty(shank_chans)
        nexttile; title(sprintf('Shank %d (Empty)', target_shank)); axis off;
        continue;
    end
    
    % 2. Sort these channels by Depth (Tip to Top)
    depths = [shank_chans.y];
    [sorted_depths, sort_idx] = sort(depths, 'ascend');
    sorted_chans = shank_chans(sort_idx);
    
    % 3. Build the Image Stack & Significance Mask
    sample_rf = sorted_chans(1).rf_map_contrast_scaleAndSmoothed;
    [~, map_width] = size(sample_rf);
    
    nChans = length(sorted_chans);
    azimuth_stack = zeros(nChans, map_width);
    sig_mask_stack = zeros(nChans, 1); 
    
    for k = 1:nChans
        map = sorted_chans(k).rf_map_contrast_scaleAndSmoothed;
        azimuth_stack(k, :) = min(map, [], 1); % Collapse Elevation (Sink)
        
        if isfield(sorted_chans(k), 'global_p_value')
            sig_mask_stack(k) = sorted_chans(k).global_p_value <= p_val_thresh;
        else
            sig_mask_stack(k) = 1; 
        end
    end
    
    % Prepare Plot Data (Invert so sink is positive)
    plot_data = azimuth_stack * -1;
    plot_data(plot_data < 0) = 0; % Clamp negative values
    
    % 4. Visualization
    nexttile;
    y_lims = [min(sorted_depths), max(sorted_depths)];
    
    % Draw the functional heatmap
    im = imagesc(x_lims, y_lims, plot_data);
    
    % Apply Alpha Mask
    alpha_vals = sig_mask_stack;
    alpha_vals(alpha_vals == 0) = nonSigAlpha;
    alpha_map = repmat(alpha_vals, 1, map_width);
    set(im, 'AlphaData', alpha_map);
    
    set(gca, 'YDir', 'normal'); 
    colormap(gca, parula);
    hold on;
    
    % =====================================================================
    % 5. ADD ANATOMICAL AXIS (The colored bar on the left)
    % =====================================================================
    % Define the physical location of the colored bar
    bar_width = 15;
    bar_left = x_lims(1) - bar_width - 5; 
    bar_right = x_lims(1) - 5;
    
    % Extract areas for this specific shank
    shank_areas = {sorted_chans.area};
    
    % Find boundaries where the area changes
    boundary_idx = [1, find(~strcmp(shank_areas(1:end-1), shank_areas(2:end))) + 1, length(shank_areas)+1];
    
    for b = 1:(length(boundary_idx)-1)
        start_idx = boundary_idx(b);
        end_idx = boundary_idx(b+1) - 1;
        
        area_name = shank_areas{start_idx};
        area_color = area_colors(area_name);
        
        % Calculate Y limits for this anatomical block
        y_start = sorted_depths(start_idx);
        y_end = sorted_depths(end_idx);
        
        % Ignore 'void' areas or very tiny slivers to keep the plot clean
        if ~strcmp(area_name, 'void')
            % Draw the colored anatomical block
            fill([bar_left, bar_right, bar_right, bar_left], ...
                 [y_start, y_start, y_end, y_end], ...
                 area_color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
             
            % Draw horizontal dashed line across the heatmap
            yline(y_start, '--', 'Color', [1 1 1 0.5], 'LineWidth', 1);
            
            % Add Text Label inside the anatomical block
            mid_y = (y_start + y_end) / 2;
            text(bar_left - 5, mid_y, area_name, 'Color', area_color, ...
                 'FontSize', 9, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end
    end
    
    % =====================================================================
    % 6. Plot RF Peaks & Gridlines
    % =====================================================================
    [~, max_col_idx] = max(plot_data, [], 2);
    deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
    peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_px;
    
    sig_indices = find(sig_mask_stack);
    if ~isempty(sig_indices)
        scatter(peak_azimuths(sig_indices), sorted_depths(sig_indices), ...
                10, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
    
    % Vertical Gridlines (Azimuth)
    xline(-120:20:120, 'w-', 'LineWidth', 1, 'Alpha', 0.3);
    
    % --- HORIZONTAL LINES AT QUARTILE DEPTHS ---
    y_span = y_lims(2) - y_lims(1);
    y_quartiles = y_lims(1) + y_span .* [0.25, 0.5, 0.75];
    yline(y_quartiles, 'r:', 'LineWidth', 1.5);
    % -------------------------------------------
    
    % Expand X-limits to fit the new anatomical bar and text!
    xlim([bar_left - 45, x_lims(2)]); 
    ax = gca;
    ax.XTick = -120:40:120;
    
    title(sprintf('Shank %d', target_shank));
    xlabel('Azimuth (deg)');
    
    if i == 1
        ylabel('Distance from tip (\mu m)');
    else
        set(gca, 'YTickLabel', []);
    end
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Sink Strength (Norm)';
sgtitle('Retinotopy by Shank with Anatomical Annotations', 'FontSize', 16, 'FontWeight', 'bold');

%% Option 2: Functional Probe Overlay (Azimuth mapped onto Anatomy)

% =========================================================================
% 1. EXTRACT PREFERRED AZIMUTH FOR EACH CHANNEL
% =========================================================================
fprintf('Extracting Preferred Azimuths...\n');
p_val_thresh = 0.05;
x_lims = [-120, 20];

% Get azimuth array based on map width
sample_rf = lfp(1).rf_map_contrast_scaleAndSmoothed;
map_width = size(sample_rf, 2);
deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
azimuths_deg = x_lims(1) + ((1:map_width) - 1) * deg_per_px;

pref_az = NaN(numel(lfp), 1);
is_sig = false(numel(lfp), 1);

for i = 1:numel(lfp)
    if isempty(lfp(i).rf_map_contrast_scaleAndSmoothed)
        continue;
    end
    
    % Check Significance
    if isfield(lfp(i), 'global_p_value')
        is_sig(i) = lfp(i).global_p_value <= p_val_thresh;
    else
        is_sig(i) = true; % Default to true if stats weren't run
    end
    
    % Find the strongest sink (minimum value) across the collapsed elevation
    map = lfp(i).rf_map_contrast_scaleAndSmoothed;
    az_profile = min(map, [], 1); 
    [~, peak_idx] = min(az_profile);
    
    pref_az(i) = azimuths_deg(peak_idx);
end

% =========================================================================
% 2. LOAD CCF VOLUMES (If not already loaded)
% =========================================================================
if ~exist('tv', 'var') || ~exist('anno', 'var') || ~exist('st', 'var')
    fprintf('Loading Allen CCF Volumes from disk...\n');
    tv = readNPY('template_volume_10um.npy'); 
    anno = readNPY('annotation_volume_10um_by_index.npy');
    st = readtable('structure_tree_safe_2017.csv');
end

% Extract Coronal Slices
ap_slice_idx = round((-insertion_AP + ccf_bregma_AP) / 10); 
coronal_slice = squeeze(tv(ap_slice_idx, :, :));
anno_slice = squeeze(anno(ap_slice_idx, :, :));

% Isolate Right Hemisphere
right_hemi_mask = false(size(anno_slice));
right_hemi_mask(:, 570:end) = true; 
anno_slice_right = anno_slice;
anno_slice_right(~right_hemi_mask) = 0; 

% =========================================================================
% 3. EXTRACT ANATOMICAL BOUNDARIES
% =========================================================================
unique_indices = unique(anno_slice_right);
unique_indices(unique_indices == 0) = []; 
target_areas = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(unique_indices)
    current_idx = unique_indices(i);
    if current_idx > height(st) || current_idx < 1, continue; end
    
    acronym = st.acronym{current_idx};
    isVisual = startsWith(acronym, 'VIS') && ~startsWith(acronym, 'VISC');
    isRSP    = startsWith(acronym, 'RSP');
    
    if isVisual || isRSP
        plot_acronym = regexprep(acronym, '[0-9/].*', ''); % Strip layer numbers
        layer_mask = (anno_slice_right == current_idx);
        
        if isKey(target_areas, plot_acronym)
            target_areas(plot_acronym) = target_areas(plot_acronym) | layer_mask;
        else
            target_areas(plot_acronym) = layer_mask;
        end
    end
end

% =========================================================================
% 4. GENERATE PLOT
% =========================================================================
fprintf('Generating Functional-Anatomical Plot...\n');
figure('Color', 'w', 'Position', [300, 100, 1000, 800]);

% A. Draw Background (Converted to True RGB to ignore the axes colormap)
% 1. Normalize the slice values between 0 and 1
slice_norm = double(coronal_slice);
slice_norm = slice_norm / max(slice_norm(:)); 
% 2. Stack 3 identical layers to create a True RGB grayscale matrix (HxWx3)
slice_rgb = cat(3, slice_norm, slice_norm, slice_norm); 

% Plot the True RGB image
imagesc(slice_rgb); hold on;

% B. Draw Faint Subcortical Boundaries
all_boundaries = edge(double(anno_slice_right), 'Sobel');
[all_y, all_x] = find(all_boundaries);
scatter(all_x, all_y, 1, [1 1 1], 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1); 

% C. Draw Bold Visual Boundaries
area_keys = keys(target_areas);
for k = 1:length(area_keys)
    area_name = area_keys{k};
    area_mask = target_areas(area_name);
    
    boundaries = edge(double(area_mask), 'Sobel');
    [edge_y, edge_x] = find(boundaries);
    scatter(edge_x, edge_y, 4, [1 1 1], 'filled', 'MarkerFaceAlpha', 0.5); % White outlines
    
    [y_coords, x_coords] = find(area_mask);
    if ~isempty(x_coords)
        text(mean(x_coords), mean(y_coords), area_name, ...
            'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Interpreter', 'none');
    end
end

% D. Plot the Functional Probe (Channels)
% Separate significant vs non-significant for cleaner plotting
sig_idx = find(is_sig & ~isnan(pref_az));
nonsig_idx = find(~is_sig & ~isnan(pref_az));

% Plot Non-Significant channels as tiny, faint grey dots
scatter([lfp(nonsig_idx).ML] / 10, [lfp(nonsig_idx).DV] / 10, ...
        8, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.3);

% Plot Significant channels colored by Preferred Azimuth
scatter([lfp(sig_idx).ML] / 10, [lfp(sig_idx).DV] / 10, ...
        30, pref_az(sig_idx), 'filled');

% E. Formatting
axis equal; 
set(gca, 'YDir', 'reverse'); 
xlim([400, 1140]); 
ylim([0, 800]);

% NOW this colormap and clim will ONLY apply to the scatter dots and colorbar!
colormap(gca, jet); 
clim(x_lims); 

cb = colorbar;
cb.Label.String = 'Preferred Azimuth (deg)';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

title(sprintf('Retinotopic Mapping across CCF Cortex (AP \\approx %d \\mu m)', insertion_AP), 'FontSize', 14);
xlabel('Medial-Lateral (CCF Pixels)'); 
ylabel('Dorsal-Ventral (CCF Pixels)');

fprintf('Done.\n');


%% Plot Azimuth as a Function of Channel Depth (8 Subplots + Anatomical Axis)

% =========================================================================
% 1. ANATOMICAL LOOKUP (Map channels to CCF Brain Areas)
% =========================================================================
fprintf('Extracting anatomical areas for each channel...\n');

% Ensure CCF Volumes are loaded
if ~exist('anno', 'var') || ~exist('st', 'var')
    error('CCF Volumes not found. Please load the Allen CCF volumes first!');
end

ap_slice_idx = round((-insertion_AP + ccf_bregma_AP) / 10); 
anno_slice = squeeze(anno(ap_slice_idx, :, :));

% Assign brain area acronyms
for i = 1:numel(lfp)
    ml_idx = round(lfp(i).ML / 10);
    dv_idx = round(lfp(i).DV / 10);
    
    if dv_idx > 0 && dv_idx <= size(anno_slice, 1) && ml_idx > 0 && ml_idx <= size(anno_slice, 2)
        pixel_id = anno_slice(dv_idx, ml_idx);
        if pixel_id > 0 && pixel_id <= height(st)
            acr = st.acronym{pixel_id};
            if startsWith(acr, 'VIS') || startsWith(acr, 'RSP')
                acr = regexprep(acr, '[0-9/].*', ''); % Strip layer numbers
            end
            lfp(i).area = acr;
        else
            lfp(i).area = 'void';
        end
    else
        lfp(i).area = 'void';
    end
end

unique_areas = unique({lfp.area});
area_colors = containers.Map(unique_areas, mat2cell(lines(length(unique_areas)), ones(length(unique_areas),1), 3));

%%
% =========================================================================
% 2. GENERATE THE 8-PANEL FUNCTIONAL PLOT
% =========================================================================
fprintf('Generating 8-Panel Functional-Anatomical Plot...\n');

% Parameters
p_val_thresh = 0.01;
nonSigAlpha = 0.3;               
x_lims = [-120, 120];         
unique_shanks = [0, 1, 2, 3]; 
nShanks = length(unique_shanks);

% Create Figure (Wider to comfortably fit 8 panels + color bars)
figure('Color', 'w', 'Position', [50, 100, 1800, 700]);
t = tiledlayout(1, nShanks * 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through each Shank
for i = 1:nShanks
    target_shank = unique_shanks(i);
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    if isempty(shank_chans)
        for empty_col = 1:2
            nexttile; title(sprintf('Shank %d\n(Empty)', target_shank)); axis off;
        end
        continue;
    end
    
    % Find unique Columns (X-coordinates) for this shank
    x_coords = round([shank_chans.x]);
    unique_cols = unique(x_coords);
    
    for col_idx = 1:length(unique_cols)
        target_col_x = unique_cols(col_idx);
        col_chans = shank_chans(x_coords == target_col_x);
        
        if isempty(col_chans)
            nexttile; title(sprintf('Shank %d | Col %d\n(Empty)', target_shank, col_idx)); axis off;
            continue;
        end
        
        % Sort channels by Depth (Tip to Top)
        depths = [col_chans.y];
        [sorted_depths, sort_idx] = sort(depths, 'ascend');
        sorted_chans = col_chans(sort_idx);
        
        % Build Image Stack
        sample_rf = sorted_chans(1).rf_map_contrast_scaleAndSmoothed;
        [~, map_width] = size(sample_rf);
        nChans = length(sorted_chans);
        azimuth_stack = zeros(nChans, map_width);
        sig_mask_stack = zeros(nChans, 1); 
        
        for k = 1:nChans
            map = sorted_chans(k).rf_map_contrast_scaleAndSmoothed;
            azimuth_stack(k, :) = min(map, [], 1); 
            
            if isfield(sorted_chans(k), 'global_p_value')
                sig_mask_stack(k) = sorted_chans(k).global_p_value <= p_val_thresh;
            else
                sig_mask_stack(k) = 1; 
            end
        end
        
        % Visualization
        nexttile;
        y_lims = [min(sorted_depths), max(sorted_depths)];
        
        plot_data = azimuth_stack * -1;
        plot_data(plot_data < 0) = 0; 
        
        im = imagesc(x_lims, y_lims, plot_data);
        
        % Alpha Mask
        alpha_vals = sig_mask_stack;
        alpha_vals(alpha_vals == 0) = nonSigAlpha;
        alpha_map = repmat(alpha_vals, 1, map_width);
        set(im, 'AlphaData', alpha_map);
        
        set(gca, 'YDir', 'normal'); 
        colormap(gca, parula);
        hold on;
        
        % =================================================================
        % 3. ADD ANATOMICAL AXIS (Colored bar on the left)
        % =================================================================
        bar_width = 15;
        bar_left = x_lims(1) - bar_width - 5; 
        bar_right = x_lims(1) - 5;
        
        col_areas = {sorted_chans.area};
        boundary_idx = [1, find(~strcmp(col_areas(1:end-1), col_areas(2:end))) + 1, length(col_areas)+1];
        
        for b = 1:(length(boundary_idx)-1)
            start_idx = boundary_idx(b);
            end_idx = boundary_idx(b+1) - 1;
            
            area_name = col_areas{start_idx};
            area_color = area_colors(area_name);
            
            y_start = sorted_depths(start_idx);
            y_end = sorted_depths(end_idx);
            
            if ~strcmp(area_name, 'void')
                % Draw the colored block
                fill([bar_left, bar_right, bar_right, bar_left], ...
                     [y_start, y_start, y_end, y_end], ...
                     area_color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                 
                % Draw horizontal dashed boundary line
                yline(y_start, '--', 'Color', [1 1 1 0.5], 'LineWidth', 1);
                
                % Add Text Label
                mid_y = (y_start + y_end) / 2;
                text(bar_left - 5, mid_y, area_name, 'Color', area_color, ...
                     'FontSize', 8, 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
            end
        end
        
        % =================================================================
        % 4. Plot RF Peaks & Formatting
        [~, max_col_idx] = max(plot_data, [], 2);
        deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
        peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_px;
        
        sig_indices = find(sig_mask_stack);
        if ~isempty(sig_indices)
            scatter(peak_azimuths(sig_indices), sorted_depths(sig_indices), ...
                    8, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        
        xline(-120:20:120, 'w-', 'LineWidth', 1, 'Alpha', 0.3);
        
        % Expand X-limits to fit the anatomical bar and labels
        xlim([bar_left - 45, x_lims(2)]); 
        ax = gca;
        ax.XTick = -120:40:120;
        
        col_label = 'Left';
        if col_idx == 2; col_label = 'Right'; end
        title(sprintf('Shank %d\n(%s Col)', target_shank, col_label), 'FontSize', 10);
        xlabel('Azimuth (deg)');
        
        if i == 1 && col_idx == 1
            ylabel('Distance from tip (\mu m)');
        else
            set(gca, 'YTickLabel', []);
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Sink Strength (Norm)';
sgtitle('High-Resolution Retinotopy (Shank Columns + Anatomy)', 'FontSize', 16, 'FontWeight', 'bold');


%%

%%% Plot Average 2D Receptive Fields (Binned by N channels, Paginated)

% =========================================================================
% 1. PARAMETERS
% =========================================================================
n_channels_to_average = 8; % Adjust bin size here
max_plots_per_fig = 40;    % Maximum number of subplots per figure

x_lims = [-120, 120];
unique_shanks = [0, 1, 2, 3];

fprintf('Generating Averaged Receptive Fields (binned by %d channels)...\n', n_channels_to_average);

% =========================================================================
% 2. BINNING & AVERAGING LOOP WITH PAGINATION (SHANK-BY-SHANK)
% =========================================================================
total_figs = 0; % Track overall generated figures

for s = 1:length(unique_shanks)
    target_shank = unique_shanks(s);
    
    % Filter channels belonging to this shank
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    if isempty(shank_chans)
        continue;
    end
    
    % Sort these channels geometrically by Depth (Top to Tip / Descending)
    % Higher y-values (shallower) come first, lowest depth (tip) comes last
    [~, sort_idx] = sort([shank_chans.y], 'descend');
    sorted_chans = shank_chans(sort_idx);
    
    % Calculate how many chunks we will get for this shank
    n_total = length(sorted_chans);
    num_chunks = ceil(n_total / n_channels_to_average);
    
    plot_idx = 0; % Keeps track of how many plots we've drawn in the current figure
    part_idx = 1; % Keeps track of the figure part number for this specific shank
    
    for c = 1:num_chunks
        
        % --- PAGINATION LOGIC: Create new figure if needed ---
        if plot_idx == 0
            figure('Color', 'w', 'Position', [50, 50, 1600, 900], 'Name', sprintf('Shank %d - Averaged RFs - Part %d', target_shank, part_idx));
            tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
            total_figs = total_figs + 1;
        end
        
        % Define the start and end indices for this chunk
        start_idx = (c-1) * n_channels_to_average + 1;
        end_idx = min(c * n_channels_to_average, n_total);
        
        chunk = sorted_chans(start_idx:end_idx);
        
        % Extract valid maps (ignoring channels that might have empty maps)
        valid_maps = {};
        for k = 1:length(chunk)
            if ~isempty(chunk(k).rf_map_contrast_scaleAndSmoothed)
                valid_maps{end+1} = chunk(k).rf_map_contrast_scaleAndSmoothed;
            end
        end
        
        if isempty(valid_maps)
            continue; % Skip if this entire chunk has no RF data
        end
        
        % 1. Calculate the Mean RF Map
        stacked_maps = cat(3, valid_maps{:});
        mean_map = mean(stacked_maps, 3);
        
        % 2. Calculate the Mean Depth
        mean_depth = mean([chunk.y]);
        
        % 3. Find the Dominant Anatomical Area (Mode)
        areas = {chunk.area};
        if ~isempty(areas)
            try
                cat_areas = categorical(areas);
                dom_area = char(mode(cat_areas));
            catch
                dom_area = areas{1}; % Fallback
            end
        else
            dom_area = 'Unknown';
        end
        
        % =================================================================
        % 3. VISUALIZATION
        % =================================================================
        nexttile;
        
        % Prepare Plot Data (Invert so sink is positive, clamp negative)
        plot_data = mean_map * -1;
        plot_data(plot_data < 0) = 0; 
        
        % Draw Image
        imagesc(x_lims, [1 size(plot_data,1)], plot_data);
        colormap(gca, parula);
        
        % Format Axes
        set(gca, 'YDir', 'normal');
        ax = gca;
        ax.FontSize = 8;
        ax.XTick = [-120, -60, 0, 60, 120];
        ax.YTick = []; % Hide Elevation Y-ticks to keep small plots clean
        
        % Dynamic Title
        title_str = sprintf('Shank %d | Depth: %.0f \\mu m\nArea: %s (n=%d)', ...
            target_shank, mean_depth, dom_area, length(valid_maps));
        
        % Color code visual/retrosplenial areas in orange
        if startsWith(dom_area, 'VIS') || startsWith(dom_area, 'RSP')
            title(title_str, 'FontSize', 9, 'Color', [0.8500 0.3250 0.0980]);
        else
            title(title_str, 'FontSize', 9, 'Color', 'k');
        end
        
        % Add faint gridline down the center (0 azimuth)
        hold on; xline(0, 'w:', 'Alpha', 0.5); hold off;
        
        % --- Increment plot counter ---
        plot_idx = plot_idx + 1;
        
        % --- If we reached the maximum plots for this figure, close it out ---
        if plot_idx >= max_plots_per_fig
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = 'Mean Sink Strength (Norm)';
            sgtitle(sprintf('Shank %d: Averaged Receptive Fields (Part %d)', target_shank, part_idx), 'FontSize', 16, 'FontWeight', 'bold');
            
            plot_idx = 0; % Reset counter for the next figure
            part_idx = part_idx + 1; % Increment part number
        end
    end
    
    % =========================================================================
    % 4. FINALIZE LAST FIGURE FOR CURRENT SHANK
    % =========================================================================
    % If the last figure for this shank didn't perfectly hit the maximum limit, 
    % we still need to add its colorbar and super-title.
    if plot_idx > 0
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = 'Mean Sink Strength (Norm)';
        sgtitle(sprintf('Shank %d: Averaged Receptive Fields (Part %d)', target_shank, part_idx), 'FontSize', 16, 'FontWeight', 'bold');
    end
end

fprintf('Done. Created %d paginated figures.\n', total_figs);