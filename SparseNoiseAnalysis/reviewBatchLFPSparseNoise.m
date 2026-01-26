%% Review batch LFP Sparse Noise

dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';
dataFiles = dir(fullfile(dataDir, 'M26001_20260122*_LFP*')); % 
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

%% plot azimuth as a function of channel depth

% Parameters
x_lims = [-120, 120]; % Visual Azimuth limits
unique_shanks = [0, 1, 2, 3]; % NP 2.0 usually has 4 shanks
nShanks = length(unique_shanks);

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1400, 600]);
t = tiledlayout(1, nShanks, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through each Shank
for i = 1:nShanks
    target_shank = unique_shanks(i);
    
    % 1. Filter channels belonging to this shank
    %    We assume lfp is a 1x384 struct array
    shank_mask = [lfp.shanks] == target_shank;
    shank_chans = lfp(shank_mask);
    
    % Check if we have data for this shank
    if isempty(shank_chans)
        nexttile; title(sprintf('Shank %d (Empty)', target_shank)); axis off;
        continue;
    end
    
    % 2. Sort these channels by Depth (Tip to Top)
    %    Extract depths only for this shank's channels
    depths = [shank_chans.y];
    [sorted_depths, sort_idx] = sort(depths, 'ascend');
    
    % 3. Build the Image Stack
    %    Get width from the first channel
    sample_rf = shank_chans(1).rf_map_contrast_scaleAndSmoothed;
    [~, map_width] = size(sample_rf);
    
    nChans = length(shank_chans);
    azimuth_stack = zeros(nChans, map_width);
    
    for k = 1:nChans
        % Get index in the shank_chans array
        orig_idx = sort_idx(k);
        
        % Get Map
        map = shank_chans(orig_idx).rf_map_contrast_scaleAndSmoothed;
        
        
        % Collapse Elevation (Min/Sink)
        azimuth_profile = min(map, [], 1);
        
        % Store
        azimuth_stack(k, :) = azimuth_profile;
    end
    
    % 4. Visualization
    nexttile;
    
    y_lims = [min(sorted_depths), max(sorted_depths)];
    
    % Multiply by -1 so Sink is Hot
    imagesc(x_lims, y_lims, azimuth_stack * -1);
    
    set(gca, 'YDir', 'normal'); 
    colormap(parula);
    title(sprintf('Shank %d', target_shank));
    xlabel('Azimuth (deg)');
    
    if i == 1
        ylabel('Depth (\mu m)');
    else
        % Remove Y-labels for inner plots to save space
        set(gca, 'YTickLabel', []);
    end
    
    % 5. Trend Line
    hold on;
    [~, max_col_idx] = max(azimuth_stack * -1, [], 2);
    
    deg_per_px = (x_lims(2) - x_lims(1)) / (map_width - 1);
    peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_px;
    
    % Smooth trend line
    smooth_peaks = movmean(peak_azimuths, 5);
    % plot(smooth_peaks, sorted_depths, 'g', 'LineWidth', 1.5);
    hold off;
end

% Global Elements
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Sink Strength (Norm)';
sgtitle('Retinotopy by Shank (Azimuth vs Depth)');