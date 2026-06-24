%% 1. Load and concatenate all units (Strictly SUA)
dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';

% 1. Grab everything that ends in _RFs.mat for this subject/session
% (This will temporarily include LFP and MUA files)
allFiles = dir(fullfile(dataDir, 'M26018_20260617*_RFs.mat')); 

% 2. Identify the files we DO NOT want
isLFP = contains({allFiles.name}, 'LFP');
isMUA = contains({allFiles.name}, 'MUA');

% 3. Filter the struct array to keep only the pure SUA files
dataFiles = allFiles(~isLFP & ~isMUA);

tic
session = struct();
for ifile = 1:numel(dataFiles)
    fprintf('Loading %s...\n', dataFiles(ifile).name);
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'units');
    session(ifile).units = units;
end
toc

units = cat(1, session.units);
fprintf('Total single units loaded: %d\n', numel(units));

%% 2. Apply Post-Hoc Soft Normalization
% Define your soft factor (e.g., 2 to 5 spikes/sec equivalent)
soft_norm_factor = 2; 

fprintf('Applying soft normalization (Factor = %d)...\n', soft_norm_factor);
for i = 1:numel(units)
    % Grab the un-normalized, smoothed RF (assuming 'smthRF' or 'contrastRF_smth')
    raw_rf = units(i).smthRF; 
    
    % Find the peak magnitude for this specific unit
    peak_val = max(abs(raw_rf(:)));
    
    % Apply soft norm and save as a new field
    units(i).soft_norm_RF = raw_rf ./ (peak_val + soft_norm_factor);
end

%% 3. Discretize Depths
allocs = cat(1, units.location);
nDepthSplits = 20;
maxVal = max(allocs(:,2));
minVal = min(allocs(:,2));
edges = linspace(minVal, maxVal, nDepthSplits+1);

% Transpose to match units array shape
depth_disc = discretize(allocs(:,2), edges)'; 

%% 4. Plot MUA Average as a Function of Shank and Depth
figure('Color', 'w', 'Position', [100, 100, 1200, 900]);
sgtitle(sprintf('MUA Soft-Normed Averaged RFs (%d Depth Bins)', nDepthSplits), 'FontSize', 16);

iplot = 0;
% Loop backwards through depths so Top of the cortex is at the Top of the figure
for idepth = nDepthSplits:-1:1
    for ishank = 0:3
        iplot = iplot + 1;
        subplot(nDepthSplits, 4, iplot);
        
        % Find all units in this physical bin
        idx = find([units.shank] == ishank & depth_disc == idepth);
        
        if ~isempty(idx)
            % Extract the newly soft-normed RFs and average them
            rf_stack = cat(3, units(idx).soft_norm_RF);
            rf_avg = nanmean(rf_stack, 3);
            
            % Plot
            imagesc(flipud(rf_avg));
            axis xy;
            colormap(parula); % Or redblue if you have custom colormaps
            
            % Optional: Lock color limits so all subplots are comparable
            % clim([0, 1]); 
            
            title(sprintf('Shank %i | Depth %i (n=%d)', ishank, idepth, length(idx)), 'FontSize', 8);
            
            % Formatting ticks: Only show labels on the outer edges to save space
            xticks(linspace(1, size(rf_avg, 2), 7));
            yticks(linspace(1, size(rf_avg, 1), 4));
            
            if idepth == 1
                xticklabels(-120:40:120);
                xlabel('Azimuth');
            else
                xticklabels([]);
            end
            
            if ishank == 0
                yticklabels(-60:40:60);
                ylabel('Elevation');
            else
                yticklabels([]);
            end
            
        else
            % If no units exist in this physical space, leave it blank
            axis off;
            title(sprintf('Shank %i | Depth %i (Empty)', ishank, idepth), 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        end
    end
end