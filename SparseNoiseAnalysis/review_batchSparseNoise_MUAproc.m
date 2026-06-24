%% Review Batch MUA Sparse Noise with Pinpoint Anatomical Mapping

% =========================================================================
% 0. LOAD BATCH MUA DATA
% =========================================================================
fprintf('Loading MUA data files...\n');

% Set your data directory here
dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis'; 
dataFiles = dir(fullfile(dataDir, '*_MUA_RFs.mat')); 

if isempty(dataFiles)
    error('No _MUA_RFs.mat files found in the specified directory. Please check the path.');
end

clear session;
tic;
for ifile = 1:numel(dataFiles)
    fprintf('  -> Loading %s...\n', dataFiles(ifile).name);
    loaded_data = load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name));
    
    % Dynamically grab the mua_sites variable (or fallback to the first variable found)
    if isfield(loaded_data, 'mua_sites')
        session(ifile).mua = loaded_data.mua_sites;
    else
        fields = fieldnames(loaded_data);
        session(ifile).mua = loaded_data.(fields{1});
    end
end
toc;

% Concatenate all loaded sessions into a single massive struct array
mua_sites = cat(2, session.mua);
fprintf('Successfully aggregated %d total MUA sites across %d sessions.\n\n', numel(mua_sites), numel(dataFiles));

% =========================================================================
% 1. ANATOMICAL LOOKUP (Pinpoint to MUA Mapping)
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
        pp_shank = floor(absolute_id / 1280);
        local_id = mod(absolute_id, 1280);
        pp_y = floor(local_id / 2) * 15;
        
        pinpoint_map(valid_count).shank = pp_shank;
        pinpoint_map(valid_count).y = pp_y;
        pinpoint_map(valid_count).area = acr;
    end
end
fprintf('  -> Reconstructed spatial map for %d Pinpoint contacts.\n', valid_count);

% --- MAP TO MUA SUBSET (Spatial Nearest-Neighbor) ---
for i = 1:numel(mua_sites)
    chan_shank = mua_sites(i).shank;          % Note: 'shank' instead of 'shanks'
    chan_y     = mua_sites(i).depth_center;   % Note: 'depth_center' instead of 'y'
    
    % Filter the Pinpoint reference map to only look at the current shank
    shank_mask = [pinpoint_map.shank] == chan_shank;
    ref_shank_data = pinpoint_map(shank_mask);
    
    if isempty(ref_shank_data)
        mua_sites(i).area = 'void';
        continue;
    end
    
    % Find the Pinpoint reference contact physically closest to this MUA depth
    [~, min_idx] = min(abs([ref_shank_data.y] - chan_y));
    
    mua_sites(i).area = ref_shank_data(min_idx).area;
end

% Generate dynamic plot colors for the unique areas detected
unique_areas = unique({mua_sites.area});
area_colors = containers.Map(unique_areas, mat2cell(lines(length(unique_areas)), ones(length(unique_areas),1), 3));

fprintf('Anatomical mapping complete.\n');


% =========================================================================
% 2. PLOT RETINOTOPY WITH ANATOMICAL AXIS (4-PANEL OVERVIEW)
% =========================================================================
fprintf('Generating Functional-Anatomical Plot...\n');

x_lims = [-120, 120];         
unique_shanks = [0, 1, 2, 3]; 
nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 700]);
tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter sites belonging to this shank
    shank_mask = [mua_sites.shank] == target_shank;
    shank_chans = mua_sites(shank_mask);
    
    if isempty(shank_chans)
        nexttile; title(sprintf('Shank %d (Empty)', target_shank)); axis off;
        continue;
    end
    
    % 2. Sort these sites by Depth (Tip to Top)
    depths = [shank_chans.depth_center];
    [sorted_depths, sort_idx] = sort(depths, 'ascend');
    sorted_chans = shank_chans(sort_idx);
    
    % 3. Build the Image Stack
    % Using 'contrastRF_smth'. Change to 'smthRF' if you prefer.
    sample_rf = sorted_chans(1).contrastRF_smth; 
    [~, map_width] = size(sample_rf);
    
    nChans = length(sorted_chans);
    azimuth_stack = zeros(nChans, map_width);
    
    for k = 1:nChans
        map = sorted_chans(k).contrastRF_smth;
        
        % MUA DIFFERENCE: Firing rates are positive peaks, so we take MAX across elevation
        azimuth_stack(k, :) = max(map, [], 1); 
    end
    
    % Normalize the stack so colors look consistent (optional but recommended for MUA)
    plot_data = azimuth_stack;
    
    % 4. Visualization
    nexttile;
    y_lims = [min(sorted_depths), max(sorted_depths)];
    
    % Draw the functional heatmap
    imagesc(x_lims, y_lims, plot_data);
    
    set(gca, 'YDir', 'normal'); 
    colormap(gca, parula);
    hold on;
    
    % =====================================================================
    % 5. ADD ANATOMICAL AXIS (The colored bar on the left)
    % =====================================================================
    bar_width = 15;
    bar_left = x_lims(1) - bar_width - 5; 
    bar_right = x_lims(1) - 5;
    
    shank_areas = {sorted_chans.area};
    boundary_idx = [1, find(~strcmp(shank_areas(1:end-1), shank_areas(2:end))) + 1, length(shank_areas)+1];
    
    for b = 1:(length(boundary_idx)-1)
        start_idx = boundary_idx(b);
        end_idx = boundary_idx(b+1) - 1;
        
        area_name = shank_areas{start_idx};
        area_color = area_colors(area_name);
        
        y_start = sorted_depths(start_idx);
        y_end = sorted_depths(end_idx);
        
        % If there's only 1 site in this area block, give it a small artificial height block
        if y_start == y_end 
            y_start = y_start - 25;
            y_end = y_end + 25;
        end
        
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
    
    scatter(peak_azimuths, sorted_depths, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Vertical Gridlines (Azimuth)
    xline(-120:20:120, 'w-', 'LineWidth', 1, 'Alpha', 0.3);
    
    % Horizontal lines at quartile depths
    y_span = y_lims(2) - y_lims(1);
    if y_span > 0
        y_quartiles = y_lims(1) + y_span .* [0.25, 0.5, 0.75];
        yline(y_quartiles, 'r:', 'LineWidth', 1.5);
    end
    
    % Expand X-limits to fit the new anatomical bar and text
    xlim([bar_left - 45, x_lims(2)]); 
    ax = gca;
    ax.XTick = -120:40:120;
    
    title(sprintf('Shank %d', target_shank));
    xlabel('Azimuth (deg)');
    
    if i == 1
        ylabel('Depth Center (\mu m)');
    else
        set(gca, 'YTickLabel', []);
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'MUA Firing Rate (Norm)';
sgtitle('MUA Retinotopy by Shank with Pinpoint Annotations', 'FontSize', 16, 'FontWeight', 'bold');

% =========================================================================
% 3. PLOT MUA AVERAGE RF AS A FUNCTION OF SHANK AND DEPTH BINS
% =========================================================================
fprintf('Generating binned MUA average RF plots...\n');

nDepthSplits = 10; % Adjust the number of physical depth bins here
all_depths = [mua_sites.depth_center];
min_depth = min(all_depths);
max_depth = max(all_depths);

% Create depth bin edges dynamically
depth_edges = linspace(min_depth, max_depth, nDepthSplits + 1);
depth_edges(1) = depth_edges(1) - 1; % Expand slightly to catch edges
depth_edges(end) = depth_edges(end) + 1;
depth_disc = discretize(all_depths, depth_edges);

figure('Color', 'w', 'Position', [100, 100, 1200, 900]);
sgtitle(sprintf('MUA Averaged Receptive Fields (%d Depth Bins)', nDepthSplits), 'FontSize', 16, 'FontWeight', 'bold');
iplot = 0;

% Loop backwards through depths so Top of the cortex is at the Top of the figure
for idepth = nDepthSplits:-1:1
    for ishank = 0:3
        iplot = iplot + 1;
        subplot(nDepthSplits, 4, iplot);
        
        % Find all MUA sites in this physical bin and shank
        idx = find([mua_sites.shank] == ishank & depth_disc == idepth);
        
        if ~isempty(idx)
            % Extract the normalized/smoothed RFs and average them
            rf_stack = cat(3, mua_sites(idx).contrastRF_smth);
            rf_avg = mean(rf_stack, 3, 'omitnan');
            
            % Plot (using flipud so elevation naturally points upwards)
            imagesc(x_lims, [-60 60], flipud(rf_avg)); 
            axis xy;
            colormap(gca, parula); 
            
            % Lock color limits if you want strict comparability
            % clim([0, max(rf_avg(:))]); 
            
            bin_start = depth_edges(idepth);
            bin_end = depth_edges(idepth+1);
            title(sprintf('Shank %i | %.0f-%.0f \\mu m (n=%d)', ishank, bin_start, bin_end, length(idx)), 'FontSize', 8);
            
            % Formatting ticks: Only show labels on outer edges to save space
            xticks(linspace(x_lims(1), x_lims(2), 7));
            yticks(linspace(-60, 60, 4));
            
            if idepth == 1
                xticklabels(-120:40:120);
                xlabel('Azimuth (deg)');
            else
                xticklabels([]);
            end
            
            if ishank == 0
                yticklabels(-60:40:60);
                ylabel('Elevation (deg)');
            else
                yticklabels([]);
            end
            
        else
            % If no sites exist in this physical space, leave it blank
            axis off;
            bin_start = depth_edges(idepth);
            bin_end = depth_edges(idepth+1);
            title(sprintf('Shank %i | %.0f-%.0f \\mu m (Empty)', ishank, bin_start, bin_end), 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        end
    end
end