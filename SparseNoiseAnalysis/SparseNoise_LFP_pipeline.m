%% Sparse noise analysis pipeline
BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';
Subject = 'M26017';
Session = '20260604';
AcquisitionsID = num2str(iacq);
SorterName = 'kilosort4';

for iacq = 0:7



    grid_size = [8, 16]; % 15 deg, normal cortex
    % grid_size = [12, 24]; % 10 deg
    % grid_size = [16, 32]; % 7.5deg, normal lgn

    %% load ephys and nidaq data info


    [ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);

    %% load LFP

    lfpDir = fullfile(ephysDir,'lfp');
    lfpMetaFile =  dir(fullfile(lfpDir,'*.meta'));
    lfpBinFile =  dir(fullfile(lfpDir,'*.bin'));

    lfpMeta = ReadMeta(fullfile(lfpMetaFile.folder, lfpMetaFile.name));

    channels = 1:384;
    targetFreq = [];
    timeRangeSamples = [];

    % less efficient, but conserves RAM
    % [lfp_data, lfp_time] = ReadBinFlexible(channels, lfpMeta, lfpBinFile.name, lfpDir, timeRangeSamples, targetFreq);
    tic
    [lfp_data_full, lfp_time_full] = ReadBinOptimized(channels, lfpMeta, lfpBinFile.name, lfpDir, [], []);
    toc
    tic
    [lfp_data, lfp_time] = ReadBinOptimized(channels, lfpMeta, lfpBinFile.name, lfpDir, [], 60);
    toc
    % for acquisitions with multiple G recordings, need to add approp. offset
    % lfp_time = lfp_time + ExpInfo-offset

    lfp = struct;
    % Get the coordinates of each channel
    [x_coords, y_coords, shanks] = GetChannelLocations(lfpMeta);

    for ichannel = 1:numel(channels)
        lfp(ichannel).data = double(lfp_data(ichannel,:));
        lfp(ichannel).data_raw = lfp_data_full(ichannel,:);
        lfp(ichannel).channel = channels(ichannel);
        lfp(ichannel).x = x_coords(channels(ichannel));
        lfp(ichannel).y = y_coords(channels(ichannel));
        lfp(ichannel).shanks = shanks(channels(ichannel));

        % --- GUARANTEE FIELD EXISTENCE ---
        lfp(ichannel).rf_map_contrast = [];
        lfp(ichannel).rf_map_contrast_scaleAndSmoothed = [];
        lfp(ichannel).rf_pvalues = [];
        lfp(ichannel).rf_mask = [];
        lfp(ichannel).is_visual_channel = NaN;
        lfp(ichannel).global_p_value = NaN;
        lfp(ichannel).rf_latency_map = [];
        lfp(ichannel).peak_latency = NaN;
    end


    %% generate stimulus matrix

    % find sparse noise bin file and minimally process
    binFileDir = fullfile(BaseDir,Subject,'ephys',Session,['Processed_',AcquisitionsID], 'bonsai_files');

    files = dir(fullfile(binFileDir, '*SparseNoise*.bin'));

    % Display found files and handle different cases
    numFiles = length(files);

    if numFiles == 0
        disp('No matching files found.');
    elseif numFiles == 1
        fprintf('Found one file: %s\n', fullfile(files(1).folder, files(1).name));
    else
        warning('Multiple matching files found:');
        for k = 1:numFiles
            fprintf('%s\n', fullfile(files(k).folder, files(k).name));
        end
        warning('Using first file.')
    end

    fileID=fopen(fullfile(files(1).folder,files(1).name)); % load the bin file

    thisBinFile=fread(fileID);
    fclose(fileID);



    % Translate stimulus into -1:1 scale
    stim_matrix = zeros(1,length(thisBinFile));
    stim_matrix(thisBinFile==0)=-1;
    stim_matrix(thisBinFile==255)=1;
    stim_matrix(thisBinFile==128)=0;
    % Make a NxM grid from the stimulus log
    stim_matrix = reshape(stim_matrix, [grid_size(1),grid_size(2), length(thisBinFile)/grid_size(1)/grid_size(2)]);
    stim_matrix = stim_matrix(:,:,1:end-1); % The last upswing should be ignored


    %% 4. Get photodiode signal for frame times (Lockout State Machine + Look-Back)
    fprintf('Step 4: Extracting stimulus frame times from NIDQ photodiode...\n');
    nidqMetaFile = dir(fullfile(nidqDir, '*.meta'));
    nidqBinFile = dir(fullfile(nidqDir, '*.bin'));
    NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
    NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name, nidqBinFile.folder);

    PD = NidqBin(2,:);
    niSampRate = str2double(NidqMeta.niSampRate);
    NidqTime = (1:size(NidqBin, 2)) ./ niSampRate;

    % Isolate photodiode signal for the sparse noise stimulus
    stimIdx = find(strcmp(ExpInfo.Stimulus, 'SparseNoise'));
    if stimIdx < height(ExpInfo)
        idx = ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx + 1);
    else
        idx = ExpInfo.smp_nidq(stimIdx):numel(NidqTime);
    end
    thisPD = PD(idx);
    thisTime = NidqTime(idx);

    % --- 1. DYNAMIC THRESHOLDS ---
    fprintf('  -> Calculating dynamic hysteresis thresholds...\n');
    pd_min = prctile(thisPD, 1);
    pd_max = prctile(thisPD, 99);
    pd_range = pd_max - pd_min;

    thresh_up = pd_min + (0.70 * pd_range);
    thresh_down = pd_min + (0.30 * pd_range);

    % --- 2. FIND ALL POTENTIAL CROSSINGS ---
    up_crossings = find(thisPD(1:end-1) < thresh_up & thisPD(2:end) >= thresh_up) + 1;
    down_crossings = find(thisPD(1:end-1) > thresh_down & thisPD(2:end) <= thresh_down) + 1;

    all_crossings = [up_crossings', ones(length(up_crossings), 1); ...
        down_crossings', -ones(length(down_crossings), 1)];
    all_crossings = sortrows(all_crossings, 1);

    idx_list = all_crossings(:, 1);
    type_list = all_crossings(:, 2);

    % --- 3. THE LOCKOUT STATE MACHINE ---
    current_state = 0;
    last_valid_time = -inf;
    refractory_period = 0.025;

    stim_times = zeros(1, length(idx_list));
    valid_count = 0;

    for i = 1:length(idx_list)
        current_idx = idx_list(i);
        t = thisTime(current_idx);
        type = type_list(i);

        if type ~= current_state && (t - last_valid_time) > refractory_period

            % --- 4. THE LOOK-BACK (Find true onset) ---
            if type == 1
                lookback_idx = find(thisPD(1:current_idx) <= thresh_down, 1, 'last');
            else
                lookback_idx = find(thisPD(1:current_idx) >= thresh_up, 1, 'last');
            end

            if ~isempty(lookback_idx)
                true_onset_time = thisTime(lookback_idx);
            else
                true_onset_time = t;
            end

            valid_count = valid_count + 1;
            stim_times(valid_count) = true_onset_time;

            current_state = type;
            last_valid_time = t;
        end
    end

    stim_times = stim_times(1:valid_count);
    fprintf('  -> Logged %d clean stimulus frames.\n', valid_count);

    %% 4.5 Strict Alignment Validation
    fprintf('Step 4.5: Validating stimulus and photodiode alignment...\n');

    nRecordedResponses = length(stim_times);
    nStimFrames = size(stim_matrix, 3);

    if nRecordedResponses ~= nStimFrames
        error('ALIGNMENT ISSUE: Detected %d photodiode pulses, but .bin file contains %d frames. Manual inspection required.', nRecordedResponses, nStimFrames);
    else
        fprintf('  -> Alignment Perfect: %d frames matched.\n', nStimFrames);
    end

    %% Define time windows and get trials
    window_sec = [-0.1, 0.3];     % Epoch window
    response_window = [0, 0.3];   % Measurement window

    % Get Time Vector & Sampling Rate (Assumes all channels share lfp_time)
    time_vec = lfp_time(:);
    dt = median(diff(time_vec));
    fs = 1/dt;
    fprintf('Detected Sampling Rate: %.2f Hz\n', fs);

    % Convert Time Windows to Samples
    samps_pre  = round(abs(window_sec(1)) * fs);
    samps_post = round(window_sec(2) * fs);
    epoch_len  = samps_pre + samps_post + 1;

    % Indices for measurement (relative to epoch start)
    resp_start_idx = samps_pre + round(response_window(1) * fs);
    resp_end_idx   = samps_pre + round(response_window(2) * fs);

    % Sanity check measurement window
    if resp_start_idx < 1, resp_start_idx = 1; end
    if resp_end_idx > epoch_len, resp_end_idx = epoch_len; end

    % Find closest LFP sample index for every stimulus time
    stim_indices = interp1(time_vec, 1:length(time_vec), stim_times, 'nearest');

    % Define valid range (so epochs don't go out of bounds)
    max_ind = length(time_vec);
    valid_mask = ~isnan(stim_indices) & ...
        stim_indices > samps_pre & ...
        stim_indices < (max_ind - samps_post);

    while numel(stim_indices) > size(stim_matrix,3)
        stim_indices(end)=[];
        valid_mask(end)=[];
    end
    valid_stim_indices = stim_indices(valid_mask);
    stim_matrix_clean  = stim_matrix(:, :, valid_mask);
    contrast_matrix    = abs(stim_matrix_clean); % 1=Flash, 0=Gray
    [nGridX, nGridY, ~] = size(contrast_matrix);

    % Cell array to store the sample indices for every grid position
    pixel_epoch_indices = cell(nGridX, nGridY);
    disp('Pre-calculating stimulus indices (Low Res)...');
    for x = 1:nGridX
        for y = 1:nGridY
            trials = find(contrast_matrix(x, y, :) == 1);
            if ~isempty(trials)
                center_idxs = valid_stim_indices(trials);
                idx_mat = center_idxs(:)' + (-samps_pre : samps_post)';
                pixel_epoch_indices{x,y} = idx_mat;
            end
        end
    end

    %% Fit RFs to every channel (Low Res Pass)
    scal_f = 10;
    sigma  = 3;
    nChannels = length(lfp);
    fprintf('Processing %d channels (Low Res fit)...\n', nChannels);
    h = waitbar(0, 'Mapping Receptive Fields...');

    for iCh = 1:nChannels
        raw_trace = lfp(iCh).data(:);

        current_rf_map = zeros(nGridX, nGridY);

        for x = 1:nGridX
            for y = 1:nGridY
                idx_mat = pixel_epoch_indices{x,y};
                if ~isempty(idx_mat)
                    epochs = raw_trace(idx_mat);
                    % Subtract the baseline in one fast matrix operation
                    epochs = epochs - mean(epochs(1:samps_pre, :), 1);
                    sta_trace = mean(epochs, 2);

                    % Store Magnitude
                    window_trace = sta_trace(resp_start_idx:resp_end_idx);
                    current_rf_map(x,y) = min(window_trace);
                end
            end
        end

        % Store Maps
        lfp(iCh).rf_map_contrast = current_rf_map;

        % Create Smoothed Magnitude Map
        rf = imresize(lfp(iCh).rf_map_contrast, scal_f);
        rf = imgaussfilt(rf, sigma);
        if min(rf(:)) ~= 0
            rf = rf ./ abs(min(rf(:)));
        end
        lfp(iCh).rf_map_contrast_scaleAndSmoothed = rf;

        if mod(iCh, 10) == 0, waitbar(iCh/nChannels, h); end
    end
    close(h);
    disp('Low Res Mapping Done.');

    %% Test significance vs circ-shifted LFP
    nShuffles = 500;
    min_shift_sec = 10;
    alpha = 0.05;
    nSamples = length(lfp(1).data);
    min_shift_samp = round(min_shift_sec * fs);

    % Pre-generate shift amounts so they are identical across channels for a given shuffle
    shift_amts = randi([min_shift_samp, nSamples - min_shift_samp], 1, nShuffles);

    fprintf('Running %d permutations...\n', nShuffles);
    h = waitbar(0, 'Shuffling Channels...');

    for iCh = 1:nChannels
        raw_trace = lfp(iCh).data(:);
        real_map = lfp(iCh).rf_map_contrast;

        % Temporary variables for parallel loop reduction
        temp_counts = zeros(nGridX, nGridY);
        temp_null_max = zeros(1, nShuffles);

        % PARFOR MAGIC: Distribute 500 shuffles across CPU cores
        parfor iShuff = 1:nShuffles
            shuffled_map = zeros(nGridX, nGridY);
            shift_amt = shift_amts(iShuff);

            for x = 1:nGridX
                for y = 1:nGridY
                    idx_mat = pixel_epoch_indices{x,y};
                    if ~isempty(idx_mat)
                        shifted_idx = mod(idx_mat + shift_amt - 1, nSamples) + 1;
                        epochs = raw_trace(idx_mat);
                        % Subtract the baseline in one fast matrix operation
                        epochs = epochs - mean(epochs(1:samps_pre, :), 1);
                        sta_trace = mean(epochs, 2);
                        shuffled_map(x,y) = min(sta_trace(resp_start_idx:resp_end_idx));
                    end
                end
            end

            temp_counts = temp_counts + (shuffled_map <= real_map);
            temp_null_max(iShuff) = min(shuffled_map(:));
        end

        % Assign back to struct
        lfp(iCh).rf_pvalues = temp_counts / nShuffles;
        lfp(iCh).rf_mask    = lfp(iCh).rf_pvalues < alpha;

        real_peak = min(real_map(:));
        n_beating_nulls = sum(temp_null_max <= real_peak);
        global_p = n_beating_nulls / nShuffles;

        lfp(iCh).is_visual_channel = global_p < alpha;
        lfp(iCh).global_p_value    = global_p;

        if mod(iCh, 10) == 0, waitbar(iCh/nChannels, h); end
    end
    close(h);
    disp('Significance testing Done.');

    %% Step 4: High-Resolution Latency Extraction
    % ---------------------------------------------------------------
    % Uses lfp_time_full and lfp(iCh).data_raw to compute precise latency
    % Only runs on Significant Channels to save time.

    if exist('lfp_time_full', 'var')
        % 1. Setup High-Res Parameters
        dt_raw = median(diff(lfp_time_full));
        fs_raw = 1/dt_raw;
        fprintf('Detected High-Res Rate: %.2f Hz\n', fs_raw);

        samps_pre_raw  = round(abs(window_sec(1)) * fs_raw);
        samps_post_raw = round(window_sec(2) * fs_raw);
        resp_start_raw = samps_pre_raw + round(response_window(1) * fs_raw);
        resp_end_raw   = samps_pre_raw + round(response_window(2) * fs_raw);

        % 2. Pre-calculate High-Res Indices
        stim_indices_raw = interp1(lfp_time_full, 1:length(lfp_time_full), stim_times, 'nearest');
        max_ind_raw = length(lfp_time_full);
        valid_mask_raw = ~isnan(stim_indices_raw) & ...
            stim_indices_raw > samps_pre_raw & ...
            stim_indices_raw < (max_ind_raw - samps_post_raw);

        % Ensure lengths match (same logic as low-res safety check)
        while numel(stim_indices_raw) > size(stim_matrix,3)
            stim_indices_raw(end)=[];
            valid_mask_raw(end)=[];
        end

        valid_stim_indices_raw = stim_indices_raw(valid_mask_raw);
        contrast_matrix_raw = abs(stim_matrix(:, :, valid_mask_raw));

        pixel_indices_raw = cell(nGridX, nGridY);
        fprintf('Pre-calculating High-Res indices...\n');
        for x = 1:nGridX
            for y = 1:nGridY
                trials = find(contrast_matrix_raw(x, y, :) == 1);
                if ~isempty(trials)
                    center_idxs = valid_stim_indices_raw(trials);
                    idx_mat = center_idxs(:)' + (-samps_pre_raw : samps_post_raw)';
                    pixel_indices_raw{x,y} = idx_mat;
                end
            end
        end

        % 3. Calculate Latency
        fprintf('Extracting precise latency for significant channels...\n');
        h = waitbar(0, 'Extracting Latency...');

        for iCh = 1:nChannels
            % Skip computation if not significant
            if lfp(iCh).is_visual_channel == 0
                lfp(iCh).peak_latency = NaN;
                lfp(iCh).rf_latency_map = nan(nGridX, nGridY); % Fill with NaNs
                continue;
            end

            raw_trace = double(lfp(iCh).data_raw(:));
            %raw_trace = lfp(iCh).data_raw(:);
            current_latency_map = nan(nGridX, nGridY);

            for x = 1:nGridX
                for y = 1:nGridY
                    idx_mat = pixel_indices_raw{x,y};
                    if ~isempty(idx_mat)
                        epochs = raw_trace(idx_mat);
                        % Subtract the baseline in one fast matrix operation
                        epochs = epochs - mean(epochs(1:samps_pre, :), 1);
                        sta_trace = mean(epochs, 2);
                        % Measure Latency in High Res Window
                        window_trace = sta_trace(resp_start_raw:resp_end_raw);
                        [~, min_idx] = min(window_trace);

                        % Convert to seconds
                        latency_samps = (resp_start_raw - samps_pre_raw - 1) + min_idx;
                        current_latency_map(x,y) = latency_samps / fs_raw;
                    end
                end
            end

            lfp(iCh).rf_latency_map = current_latency_map;

            % Extract Peak Latency from the "hottest" pixel in the (Low Res) Magnitude map
            [~, best_flat_idx] = min(lfp(iCh).rf_map_contrast(:));
            lfp(iCh).peak_latency = current_latency_map(best_flat_idx);

            if mod(iCh, 10) == 0, waitbar(iCh/nChannels, h); end
        end
        close(h);
        disp('High-Res Latency Extraction Done.');
    else
        warning('lfp_time_full not found. Skipping high-res latency.');
    end


    %% save data

    lfp_save = rmfield(lfp, {'data', 'data_raw'});

    fname = [Subject, '_', Session, '_', AcquisitionsID, '_LFP_RFs_wLatency.mat'];
    try
        save(fname, 'lfp_save'); % Saves instantly!
    catch
        save(fname, 'lfp_save', '-v7.3')
    end


end

%% plot lfp for all channels, sorted by depth
%
% [depth, depth_idx] = sort([lfp.y]);
% depth_idx = fliplr(depth_idx);
% scal_f=10;
% sigma=3;
% %
% doplot=true;
% if doplot
%     % idx=1:numel(units);
%     tic
%     for ichannel = 1:384
%         % For every new set of panels, open a new figure (except for the first set)
%         if mod(ichannel-1, 35) == 0
%             figure;         % Open a new figure
%         end
%
%         % Create subplot in a 4x5 grid
%         subplot(7, 5, mod(ichannel-1, 35) + 1)
%
%         % get the receptive field
%         this_chan = depth_idx(ichannel);
%         rf = imresize(lfp(this_chan).rf_map_contrast,scal_f);
%         rf = imgaussfilt(rf,sigma);
%         rf = rf./abs(max(rf(:)));
%
%         % plot the receptive field
%         imagesc(flipud(rf))
%         colormap(gray)
%         xlabel('Azimuth')
%         ylabel('Elevation')
%         xticks(linspace(1, size(rf, 2), 13))
%         xticklabels(-120:20:120)
%         yticks(linspace(1, size(rf, 1), 7))
%         yticklabels(-60:20:60)
%         title(sprintf('Depth %i', lfp(this_chan).y))
%         colorbar
%         set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 8)
%         % axis equal,
%         axis xy
%         % colormap(redblue), clim([-abv abv])
%     end
%
%     toc
% end
%

%% plot azimuth as a function of channel depth
% % Setup Parameters
% scal_f = 10;      % Resize factor
% sigma  = 3;       % Smoothing sigma
% normalize = true; % Normalize to peak sink
%
% % Get output dimensions for resizing
% dummy_map = imresize(lfp(1).rf_map_contrast, scal_f);
% [~, new_w] = size(dummy_map);
%
% % Extract Depths and Sort
% all_depths = [lfp.y];
% [sorted_depths, sort_idx] = sort(all_depths, 'ascend');
%
% % Pre-allocate Matrices
% nChannels = length(lfp);
% azimuth_depth_sorted = zeros(nChannels, new_w);
% significance_mask_sorted = zeros(nChannels, 1); % To track which rows are sig
%
% % Process Channels in Sorted Order
% fprintf('Smoothing and sorting %d channels...\n', nChannels);
%
% for k = 1:nChannels
%     % Get the original channel index corresponding to this depth rank
%     orig_idx = sort_idx(k);
%
%     % Get Raw RF
%     raw_rf = lfp(orig_idx).rf_map_contrast;
%
%     % Resize & Smooth
%     rf_resized = imresize(raw_rf, scal_f);
%     rf_smooth = imgaussfilt(rf_resized, sigma);
%
%     % Normalize (to min/sink)
%     if normalize
%         peak_val = min(rf_smooth(:));
%         if peak_val ~= 0
%             rf_smooth = rf_smooth ./ abs(peak_val);
%         end
%     end
%
%     % Collapse Elevation to get Azimuth Profile
%     azimuth_profile = min(rf_smooth, [], 1);
%
%     % Store data
%     azimuth_depth_sorted(k, :) = azimuth_profile;
%
%     % Store significance status (computed in previous permutation step)
%     % If the field doesn't exist, assume all are sig (safety check)
%     if isfield(lfp, 'is_visual_channel')
%         significance_mask_sorted(k) = lfp(orig_idx).is_visual_channel;
%     else
%         significance_mask_sorted(k) = 1;
%     end
% end
%
% % Prepare Data for Plotting
% plot_data = azimuth_depth_sorted * -1; % Invert so sink is positive
% % Clamp negatives (sources/rebounds) to 0 so the colormap starts at 0
% plot_data(plot_data < 0) = 0;
%
% % Visualization
% figure('Color', 'w'); % White background
%
% % Define Bounds
% x_lims = [-120, 120];
% y_lims = [min(sorted_depths), max(sorted_depths)];
%
% % Create Image Object
% im = imagesc(x_lims, y_lims, plot_data);
%
% % --- KEY STEP: Handle Significance (White for Non-Sig) ---
% % We use AlphaData to make non-significant rows transparent.
% % Since the figure background is white, transparent = white.
%
% % 1. Create a full-size alpha map (nChannels x nWidth)
% % Replicate the column vector across the width
% alpha_map = repmat(significance_mask_sorted, 1, new_w);
%
% % 2. Apply Alpha Data
% % 1 = Opaque (Significant), 0 = Transparent (White/Non-Sig)
% set(im, 'AlphaData', alpha_map);
%
% % Formatting
% set(gca, 'YDir', 'normal');
% colormap(parula); % Changed to Parula
% c = colorbar;
% c.Label.String = 'Normalized Response (Sink Strength)';
%
% xlabel('Azimuth (Degrees)');
% ylabel('Probe Depth (\mu m)');
% title('Retinotopy: Azimuth vs. Physical Depth (Significant Channels Only)');
%
% % Optional: Track Peak (Only for significant channels)
% hold on;
% [~, max_col_idx] = max(plot_data, [], 2);
% deg_per_pixel = (x_lims(2) - x_lims(1)) / (new_w - 1);
% peak_azimuths = x_lims(1) + (max_col_idx - 1) * deg_per_pixel;
%
% % Filter peak line to show only significant points
% sig_indices = find(significance_mask_sorted);
% % Use dots ('.') instead of a line so gaps appear where channels are non-sig
% plot(peak_azimuths(sig_indices), sorted_depths(sig_indices), 'k.', 'MarkerSize', 8);
% hold off;