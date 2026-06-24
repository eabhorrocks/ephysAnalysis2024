%% Sparse noise analysis pipeline
BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';
Subject = 'M26018';
Session = '20260621';
SorterName = 'kilosort4';

% grid_size = [8, 16]; % 15 deg, normal cortex
% grid_size = [12, 24]; % 10 deg
% grid_size = [16, 32]; % 7.5deg, normal lgn
grid_size = [12, 30]; % SS19 rig

% nidq photodiode processing parameters
refractory_period = 0.00;
down_thresh = 0.5;
up_thresh = 0.5;
upper_centile = 95;
lower_centile = 5;
maxJit = 1.75;

override_photodiode = true;     % Set to true to bypass PD extraction
expected_frame_dur  = 0.1;      % Fixed duration between frames (in seconds)

for iacq = 0:5
    AcquisitionsID = num2str(iacq);

    %% load ephys and nidaq data info
    [ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session, AcquisitionsID);

    %% load LFP
    lfpDir = fullfile(ephysDir, 'lfp');
    lfpMetaFile =  dir(fullfile(lfpDir, '*.meta'));
    lfpBinFile =  dir(fullfile(lfpDir, '*.bin'));
    lfpMeta = ReadMeta(fullfile(lfpMetaFile.folder, lfpMetaFile.name));
    channels = 1:384;

    tic
    [lfp_data_full, lfp_time_full] = ReadBinOptimized(channels, lfpMeta, lfpBinFile.name, lfpDir, [], []);
    toc
    tic
    [lfp_data, lfp_time] = ReadBinOptimized(channels, lfpMeta, lfpBinFile.name, lfpDir, [], 60);
    toc

    lfp = struct;
    % Get the coordinates of each channel
    [x_coords, y_coords, shanks] = GetChannelLocations(lfpMeta);

    for ichannel = 1:numel(channels)
        lfp(ichannel).data = double(lfp_data(ichannel, :));
        % lfp_data_full is explicitly omitted here to save RAM. Read directly later.

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
    binFileDir = fullfile(BaseDir, Subject, 'ephys', Session, ['Processed_', AcquisitionsID], 'bonsai_files');
    files = dir(fullfile(binFileDir, '*SparseNoise*.bin'));

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

    fileID = fopen(fullfile(files(1).folder, files(1).name));
    thisBinFile = fread(fileID);
    fclose(fileID);

    % Translate stimulus into -1:1 scale
    stim_matrix = zeros(1, length(thisBinFile));
    stim_matrix(thisBinFile == 0) = -1;
    stim_matrix(thisBinFile == 255) = 1;
    stim_matrix(thisBinFile == 128) = 0;

    % Make a NxM grid from the stimulus log
    stim_matrix = reshape(stim_matrix, [grid_size(1), grid_size(2), length(thisBinFile)/grid_size(1)/grid_size(2)]);
    stim_matrix = stim_matrix(:, :, 1:end-1);

    %% 4. Get photodiode signal for frame times
    fprintf('Step 4: Extracting stimulus frame times from NIDQ photodiode...\n');
    

    
    nidqMetaFile = dir(fullfile(nidqDir, '*.meta'));
    nidqBinFile = dir(fullfile(nidqDir, '*.bin'));
    NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
    NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name, nidqBinFile.folder);
    PD = NidqBin(2, :);
    niSampRate = str2double(NidqMeta.niSampRate);
    NidqTime = (1:size(NidqBin, 2)) ./ niSampRate;
    
    stimIdx = find(strcmp(ExpInfo.Stimulus, 'SparseNoise'));
    if stimIdx < height(ExpInfo)
        idx = ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx + 1);
    else
        idx = ExpInfo.smp_nidq(stimIdx):numel(NidqTime);
    end
    thisPD = PD(idx);
    thisTime = NidqTime(idx);
    
    nStimFrames = size(stim_matrix, 3);
    
    if override_photodiode
        fprintf('  -> OVERRIDE ENABLED: Assuming fixed %.3fs spacing.\n', expected_frame_dur);
        
        % We still need the VERY FIRST pulse to anchor the timeline
        pd_min = prctile(thisPD, lower_centile);
        pd_max = prctile(thisPD, upper_centile);
        thresh_up = pd_min + (up_thresh * (pd_max - pd_min));
        
        first_up_idx = find(thisPD(1:end-1) < thresh_up & thisPD(2:end) >= thresh_up, 1, 'first');
        
        if isempty(first_up_idx)
            warning('Could not find even a single PD pulse! Anchoring to the start of the epoch.');
            anchor_time = thisTime(1);
        else
            anchor_time = thisTime(first_up_idx + 1);
            fprintf('  -> Anchored timeline to first PD pulse at t = %.3fs\n', anchor_time);
        end
        
        % Mathematically generate the perfect timeline
        stim_times = anchor_time + (0 : nStimFrames - 1) * expected_frame_dur;
        valid_count = length(stim_times);
        
    else
        % --- ORIGINAL PHOTODIODE EXTRACTION LOGIC ---
        fprintf('  -> Calculating dynamic hysteresis thresholds...\n');
        pd_min = prctile(thisPD, lower_centile);
        pd_max = prctile(thisPD, upper_centile);
        pd_range = pd_max - pd_min;
        thresh_up = pd_min + (up_thresh * pd_range);
        thresh_down = pd_min + (down_thresh * pd_range);
        
        up_crossings = find(thisPD(1:end-1) < thresh_up & thisPD(2:end) >= thresh_up) + 1;
        down_crossings = find(thisPD(1:end-1) > thresh_down & thisPD(2:end) <= thresh_down) + 1;
        all_crossings = [up_crossings', ones(length(up_crossings), 1); ...
            down_crossings', -ones(length(down_crossings), 1)];
        all_crossings = sortrows(all_crossings, 1);
        idx_list = all_crossings(:, 1);
        type_list = all_crossings(:, 2);
        
        current_state = 0;
        last_valid_time = -inf;
        stim_times = zeros(1, length(idx_list));
        valid_count = 0;
        
        for i = 1:length(idx_list)
            current_idx = idx_list(i);
            t = thisTime(current_idx);
            type = type_list(i);
            if type ~= current_state && (t - last_valid_time) > refractory_period
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
    end
    
    fprintf('  -> Logged %d clean stimulus frames.\n', valid_count);

    %% Auto-repair block for potentially completely missed frames

    stim_times = stim_times(1:valid_count);
    fprintf('  -> Logged %d clean stimulus frames.\n', valid_count);

    %% 4.2 Automated Monitor Freeze Repair
    % Detect if frames are missing, and if it's strictly due to a frozen monitor
    nStimFrames = size(stim_matrix, 3);
    nRecordedResponses = length(stim_times);

    if nRecordedResponses < nStimFrames
        fprintf('  -> WARNING: Missing %d frames. Checking if monitor froze...\n', nStimFrames - nRecordedResponses);
        
        expected_dt = median(diff(stim_times));
        dt = diff(stim_times);
        
        % Ignore standard monitor jitter (up to 1.75x) and only flag actual dropped frames
        dropped_idx = find(dt > maxJit * expected_dt);        
        is_clean_freeze = true;
        total_missing_calculated = 0;
        
        for k = 1:length(dropped_idx)
            gap_ratio = dt(dropped_idx(k)) / expected_dt;
            
            % Check if the gap is very close to a clean integer (e.g. 2.0x, 3.0x expected duration)
            if abs(gap_ratio - round(gap_ratio)) < 0.15 
                total_missing_calculated = total_missing_calculated + (round(gap_ratio) - 1);
            else
                is_clean_freeze = false; % Gap is messy, likely a thresholding issue
            end
        end
        
        % ONLY interpolate if every gap is a clean integer AND the math perfectly solves the missing frames
        if is_clean_freeze && (total_missing_calculated == (nStimFrames - nRecordedResponses))
            fprintf('  -> Confirmed exactly %d dropped monitor frames. Safely applying interpolation...\n', total_missing_calculated);
            
            repaired_stim_times = [];
            last_idx = 1;
            for k = 1:length(dropped_idx)
                gap_idx = dropped_idx(k);
                repaired_stim_times = [repaired_stim_times, stim_times(last_idx:gap_idx)];
                
                n_missing = round(dt(gap_idx) / expected_dt) - 1; 
                t_start = stim_times(gap_idx);
                t_end = stim_times(gap_idx + 1);
                
                inserted_times = linspace(t_start, t_end, n_missing + 2);
                repaired_stim_times = [repaired_stim_times, inserted_times(2:end-1)];
                
                last_idx = gap_idx + 1;
            end
            repaired_stim_times = [repaired_stim_times, stim_times(last_idx:end)];
            stim_times = repaired_stim_times;
            
            fprintf('  -> SUCCESS: Timeline mathematically repaired. Total frames now match (%d).\n', length(stim_times));
        else
            fprintf('  -> REPAIR ABORTED: Gaps are irregular or math does not align. This is a Thresholding Issue.\n');
            % Do nothing. Let Step 4.5 throw the hard error so you can tune your thresholds.
        end
    elseif nRecordedResponses > nStimFrames
        fprintf('  -> WARNING: Detected MORE frames than expected. This is a Thresholding Issue.\n');
    end

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
    window_sec = [-0.1, 0.3];
    response_window = [0, 0.3];
    time_vec = lfp_time(:);
    dt = median(diff(time_vec));
    fs = 1/dt;
    fprintf('Detected Sampling Rate: %.2f Hz\n', fs);

    samps_pre  = round(abs(window_sec(1)) * fs);
    samps_post = round(window_sec(2) * fs);
    epoch_len  = samps_pre + samps_post + 1;

    % FIX: +1 for accurate MATLAB index alignment
    resp_start_idx = samps_pre + 1 + round(response_window(1) * fs);
    resp_end_idx   = samps_pre + 1 + round(response_window(2) * fs);

    if resp_start_idx < 1, resp_start_idx = 1; end
    if resp_end_idx > epoch_len, resp_end_idx = epoch_len; end

    stim_indices = interp1(time_vec, 1:length(time_vec), stim_times, 'nearest');
    max_ind = length(time_vec);
    valid_mask = ~isnan(stim_indices) & ...
        stim_indices > samps_pre & ...
        stim_indices < (max_ind - samps_post);

    while numel(stim_indices) > size(stim_matrix, 3)
        stim_indices(end) = [];
        valid_mask(end) = [];
    end

    valid_stim_indices = stim_indices(valid_mask);
    stim_matrix_clean  = stim_matrix(:, :, valid_mask);
    contrast_matrix    = abs(stim_matrix_clean);
    [nGridX, nGridY, ~] = size(contrast_matrix);
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
                    epochs = epochs - mean(epochs(1:samps_pre, :), 1);
                    sta_trace = mean(epochs, 2);
                    window_trace = sta_trace(resp_start_idx:resp_end_idx);
                    current_rf_map(x,y) = min(window_trace);
                end
            end
        end
        lfp(iCh).rf_map_contrast = current_rf_map;
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
    nShuffles = 1;
    min_shift_sec = 10;
    alpha = 0.05;
    nSamples = length(lfp(1).data);
    min_shift_samp = round(min_shift_sec * fs);
    shift_amts = randi([min_shift_samp, nSamples - min_shift_samp], 1, nShuffles);

    fprintf('Running %d permutations across %d channels...\n', nShuffles, nChannels);

    % Pre-allocate outputs to avoid parfor structure slicing issues
    out_pvalues  = cell(nChannels, 1);
    out_mask     = cell(nChannels, 1);
    out_is_vis   = cell(nChannels, 1);
    out_global_p = cell(nChannels, 1);

    % FIX: Parallelize over channels (spins up the pool exactly once)
    parfor iCh = 1:nChannels
        raw_trace = lfp(iCh).data(:);
        real_map = lfp(iCh).rf_map_contrast;
        real_peak = min(real_map(:));

        temp_counts = zeros(nGridX, nGridY);
        temp_null_max = zeros(1, nShuffles);

        for iShuff = 1:nShuffles
            shuffled_map = zeros(nGridX, nGridY);
            shift_amt = shift_amts(iShuff);

            for x = 1:nGridX
                for y = 1:nGridY
                    idx_mat = pixel_epoch_indices{x,y};
                    if ~isempty(idx_mat)
                        shifted_idx = mod(idx_mat + shift_amt - 1, nSamples) + 1;
                        epochs = raw_trace(shifted_idx); % Correct shifted index
                        epochs = epochs - mean(epochs(1:samps_pre, :), 1);
                        sta_trace = mean(epochs, 2);
                        shuffled_map(x,y) = min(sta_trace(resp_start_idx:resp_end_idx));
                    end
                end
            end
            temp_counts = temp_counts + (shuffled_map <= real_map);
            temp_null_max(iShuff) = min(shuffled_map(:));
        end

        out_pvalues{iCh} = temp_counts / nShuffles;
        out_mask{iCh}    = (temp_counts / nShuffles) < alpha;

        n_beating_nulls = sum(temp_null_max <= real_peak);
        global_p = n_beating_nulls / nShuffles;

        out_is_vis{iCh}   = global_p < alpha;
        out_global_p{iCh} = global_p;
    end

    for iCh = 1:nChannels
        lfp(iCh).rf_pvalues = out_pvalues{iCh};
        lfp(iCh).rf_mask    = out_mask{iCh};
        lfp(iCh).is_visual_channel = out_is_vis{iCh};
        lfp(iCh).global_p_value    = out_global_p{iCh};
    end
    disp('Significance testing Done.');

    %% Step 4: High-Resolution Latency Extraction
    if exist('lfp_time_full', 'var')
        dt_raw = median(diff(lfp_time_full));
        fs_raw = 1/dt_raw;
        fprintf('Detected High-Res Rate: %.2f Hz\n', fs_raw);

        samps_pre_raw  = round(abs(window_sec(1)) * fs_raw);
        samps_post_raw = round(window_sec(2) * fs_raw);

        % FIX: +1 for accurate MATLAB index alignment
        resp_start_raw = samps_pre_raw + 1 + round(response_window(1) * fs_raw);
        resp_end_raw   = samps_pre_raw + 1 + round(response_window(2) * fs_raw);

        stim_indices_raw = interp1(lfp_time_full, 1:length(lfp_time_full), stim_times, 'nearest');
        max_ind_raw = length(lfp_time_full);
        valid_mask_raw = ~isnan(stim_indices_raw) & ...
            stim_indices_raw > samps_pre_raw & ...
            stim_indices_raw < (max_ind_raw - samps_post_raw);

        while numel(stim_indices_raw) > size(stim_matrix, 3)
            stim_indices_raw(end) = [];
            valid_mask_raw(end) = [];
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

        fprintf('Extracting precise latency for significant channels...\n');
        h = waitbar(0, 'Extracting Latency...');

        for iCh = 1:nChannels
            if lfp(iCh).is_visual_channel == 0
                lfp(iCh).peak_latency = NaN;
                lfp(iCh).rf_latency_map = nan(nGridX, nGridY);
                continue;
            end

            % FIX: Pull direct from matrix memory, not struct
            raw_trace = double(lfp_data_full(iCh, :));
            raw_trace = raw_trace(:);

            current_latency_map = nan(nGridX, nGridY);
            for x = 1:nGridX
                for y = 1:nGridY
                    idx_mat = pixel_indices_raw{x,y};
                    if ~isempty(idx_mat)
                        epochs = raw_trace(idx_mat);

                        % FIX: Ensure using high-res pre-samples for baseline
                        epochs = epochs - mean(epochs(1:samps_pre_raw, :), 1);
                        sta_trace = mean(epochs, 2);

                        window_trace = sta_trace(resp_start_raw:resp_end_raw);
                        [~, min_idx] = min(window_trace);

                        % FIX: Accurate latency math
                        idx_within_epoch = (resp_start_raw - 1) + min_idx;
                        stimulus_idx = samps_pre_raw + 1;
                        latency_samps = idx_within_epoch - stimulus_idx;

                        current_latency_map(x,y) = latency_samps / fs_raw;
                    end
                end
            end

            lfp(iCh).rf_latency_map = current_latency_map;
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
    % data_raw is gone, so just clear 'data'
    lfp_save = rmfield(lfp, 'data');
    fname = [Subject, '_', Session, '_', AcquisitionsID, '_LFP_RFs_wLatency.mat'];
    try
        save(fname, 'lfp_save');
    catch
        save(fname, 'lfp_save', '-v7.3')
    end
end


%% debugging

doDebug=false
if doDebug
% 1. Calculate inter-frame intervals
dt = diff(stim_times);

% Plot Histogram
figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
subplot(1, 2, 1);
histogram(dt, 100, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Time between frames (s)');
ylabel('Count');
title('Histogram of Inter-Frame Intervals');
grid on;

% 2. Find the maximum gap
[max_dt, max_idx] = max(dt);
t1 = stim_times(max_idx);
t2 = stim_times(max_idx + 1);

fprintf('Maximum inter-frame interval is %.4f seconds.\n', max_dt);
fprintf('Occurred between Frame %d (t=%.2fs) and Frame %d (t=%.2fs).\n', ...
    max_idx, t1, max_idx + 1, t2);

% 3. Plot the raw photodiode signal around this gap
% Define a window to look at (e.g., 0.2 seconds before the gap and 0.2s after)
window_pad = 0.2; 
window_start = t1 - window_pad;
window_end = t2 + window_pad;

% Find the indices in your raw NIDQ time vector that fall in this window
plot_idx = thisTime >= window_start & thisTime <= window_end;

subplot(1, 2, 2);
plot(thisTime(plot_idx), thisPD(plot_idx), 'k', 'LineWidth', 1.5);
hold on;

% Draw vertical red lines where the code ACTUALLY detected the frames
xline(t1, 'r--', sprintf('Frame %d', max_idx), 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Color', 'r', 'LineWidth', 1.5);
xline(t2, 'r--', sprintf('Frame %d', max_idx+1), 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'Color', 'r', 'LineWidth', 1.5);

% If you still have your thresholds in the workspace, plot them too!
if exist('thresh_up', 'var') && exist('thresh_down', 'var')
    yline(thresh_up, 'b:', 'Thresh Up', 'LineWidth', 1);
    yline(thresh_down, 'b:', 'Thresh Down', 'LineWidth', 1);
end

xlabel('NIDQ Time (s)');
ylabel('Photodiode Voltage (V)');
title(sprintf('Raw Signal Around Max Gap (%.3f s)', max_dt));
grid on;
end