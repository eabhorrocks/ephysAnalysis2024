%% analyse receptive fields from multiple acquisitons

%% load processed acquisition .mat files and make allUnits var
dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';
searchPattern = 'M26001_*_RFs.mat'; % Matches all M26001 RF files
fileList = dir(fullfile(dataDir, searchPattern));
allUnits = [];
fprintf('Found %d files. Starting loading...\n', length(fileList));

for i = 1:length(fileList)
    currentFileName = fileList(i).name;
    fullPath = fullfile(dataDir, currentFileName);
    load(fullPath, 'units');
    allUnits = cat(1,allUnits, units);
end

% 5. Summary
fprintf('------------------------------------------\n');
fprintf('Total units loaded: %d\n', length(allUnits));


%% cycle through all receptiove fields and label

% --- SETTINGS ---
fr_thresh = 2.0; % Set your threshold here (e.g., 1 Hz)

% 1. Initialize field
if ~isfield(allUnits, 'isGoodRF')
    [allUnits.isGoodRF] = deal(nan);
end

% 2. Identify units to curate vs skip
% Assuming 'firing_rate' exists. If not, create dummy data to prevent crash.
if ~isfield(allUnits, 'firing_rate')
    error('Field "firing_rate" is missing from allUnits.');
end

% Get list of indices that pass the threshold
validIndices = find([allUnits.firing_rate] >= fr_thresh);
skippedIndices = find([allUnits.firing_rate] < fr_thresh);

% Automatically mark skipped units as BAD (0)
for s = skippedIndices
    allUnits(s).isGoodRF = 0;
end
fprintf('Auto-marked %d units as Bad (FR < %.1f Hz).\n', length(skippedIndices), fr_thresh);

% 3. Create Figure
hFig = figure('Name', 'RF Curator', 'NumberTitle', 'off');
fprintf('Starting curation on %d remaining units...\n', length(validIndices));
fprintf('Controls: [y] = Good, [n] = Bad, [q] = Quit, [backspace] = Go back\n');

% 4. Loop using an index pointer 'k' for the validIndices list
k = 180; 

while k <= length(validIndices)
    
    % Get the actual unit index from our filtered list
    iunit = validIndices(k);
    
    % --- PLOTTING ---
    figure(hFig); clf;
    set(hFig, 'CurrentCharacter', ' '); % Reset buffer
    
    thisMap_s = allUnits(iunit).normRF;
    frStr = sprintf('%.1f Hz', allUnits(iunit).firing_rate);

    if isempty(thisMap_s)
        text(0.5, 0.5, 'Empty RF Data', 'HorizontalAlignment', 'center');
    else
        imagesc(flipud(thisMap_s));
        xlabel('Azimuth');
        ylabel('Elevation');
        xticks(linspace(1, size(thisMap_s, 2), 13));
        xticklabels(-120:20:120);
        yticks(linspace(1, size(thisMap_s, 1), 7));
        yticklabels(-60:20:60);
        colorbar;
        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 8);
        axis xy
    end
    
    % Status
    statusStr = 'Pending';
    if ~isnan(allUnits(iunit).isGoodRF)
        if allUnits(iunit).isGoodRF == 1, statusStr = 'GOOD (y)'; end
        if allUnits(iunit).isGoodRF == 0, statusStr = 'BAD (n)'; end
    end
    
    title(sprintf('Unit %d (Index %d/%d) | FR: %s | Status: %s\nPress [y] or [n]', ...
        iunit, k, length(validIndices), frStr, statusStr), 'Interpreter', 'none');

    % --- INTERACTION ---
    validInput = false;
    while ~validInput
        try
            wasKeyPress = waitforbuttonpress;
            if wasKeyPress == 1
                key = get(hFig, 'CurrentCharacter');
                switch lower(key)
                    case 'y'
                        allUnits(iunit).isGoodRF = 1;
                        fprintf('Unit %d: Good (y)\n', iunit);
                        k = k + 1; % Move to next valid unit
                        validInput = true;
                    case 'n'
                        allUnits(iunit).isGoodRF = 0;
                        fprintf('Unit %d: Bad (n)\n', iunit);
                        k = k + 1; % Move to next valid unit
                        validInput = true;
                    case 8 % Backspace
                        if k > 1
                            k = k - 1; % Move to previous valid unit
                            fprintf('Back to Unit %d\n', validIndices(k));
                        end
                        validInput = true;
                    case 'q'
                        close(hFig);
                        fprintf('Quitting.\n');
                        return;
                end
            end
        catch
            return;
        end
    end
end

close(hFig);
fprintf('Done! All passing units curated.\n');

%% get rf locations

% 1. Define your visual boundaries (matches your plotting code)
az_lims = [-120 120]; % Azimuth limits
el_lims = [-60 60];   % Elevation limits

% 2. Loop through all units
for i = 1:length(allUnits)
    
    rf = allUnits(i).normRF;
    
    if ~isempty(rf) && allUnits(i).isGoodRF == 1 % Only process existing (and optionally GOOD) RFs
        
        % A. Find the linear index of the maximum value
        [maxVal, maxLinearIdx] = max(rf, [], 'all', 'linear');
        
        % B. Convert linear index to subscripts (row, col)
        [r, c] = ind2sub(size(rf), maxLinearIdx);
        
        % C. Store raw matrix indices
        allUnits(i).rf_max_row = r;
        allUnits(i).rf_max_col = c;
        
        % D. Convert to Degrees (Space)
        % Create vector of positions for Azimuth (columns)
        az_vec = linspace(az_lims(1), az_lims(2), size(rf, 2));
        
        % Create vector of positions for Elevation (rows)
        % Note: Since you use flipud() to plot, Row 1 is the Bottom (-60)
        el_vec = linspace(el_lims(1), el_lims(2), size(rf, 1));
        
        % Map indices to degrees
        allUnits(i).rf_azimuth = az_vec(c);
        allUnits(i).rf_elevation = el_vec(r);
        
    else
        % Fill with NaNs if no RF exists or it was marked bad
        allUnits(i).rf_max_row = nan;
        allUnits(i).rf_max_col = nan;
        allUnits(i).rf_azimuth = nan;
        allUnits(i).rf_elevation = nan;
    end
end

fprintf('Done. Added rf_azimuth and rf_elevation to allUnits.\n');


%% plot unit locations and mark good rfs

good_idx = find([allUnits.isGoodRF]);

figure, hold on
plot(allLocs(:,1), allLocs(:,2), 'k.')
plot(allLocs(good_idx,1), allLocs(good_idx,2), 'r*')

%% 

sh3_units_idx = allLocs(:,1)>200 ;

sh3Units = allUnits(sh3_units_idx & [allUnits.isGoodRF]');
sh3depths = cat(1,sh3Units.location);
sh3depths = sh3depths(:,2);

plot([sh3Units.rf_azimuth], sh3depths, 'k.')

%%

allDepths = allLocs(:,2);

plot([allUnits.rf_azimuth], allDepths, 'k.')
