% Batch processing script for receptive field mapping (MUA / Spatial Grouping)
% =========================================================================
% ==                  DEFINE YOUR RECORDINGS HERE                        ==
% =========================================================================
% Add each unique recording as an element in the 'recordings' struct array.
recordings = struct([]); % Initialize the struct array
recordings(1).Subject = 'M26018';
recordings(1).Session = '20260617';
recordings(1).AcquisitionsID = '0';
recordings(1).SorterName = 'kilosort4';
for irecording = 2:8
    recordings(irecording) = recordings(1);
    recordings(irecording).AcquisitionsID = [num2str(irecording-1)];
end

% =========================================================================
% Main directory for all subject data
BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';

% RF mapping parameters
grid_size = [12, 24]; % 15 deg, normal cortex
scal_f = 10; % Scaling factor for RF upsampling
sigma = 3;   % Sigma for Gaussian smoothing of the RF

% MUA Spatial Grouping Parameters
mua_depth_bin_size = 50; % Spatial bin size for depth discretization (micrometers)

% =========================================================================
% ==                      START BATCH PROCESSING LOOP                    ==
% =========================================================================
for iRec = 3:length(recordings)
    
    % --- Get parameters for the current recording ---
    Subject = recordings(iRec).Subject;
    Session = recordings(iRec).Session;
    AcquisitionsID = recordings(iRec).AcquisitionsID;
    SorterName = recordings(iRec).SorterName;
    
    fprintf('\n\n============================================================\n');
    fprintf('Processing recording %d of %d (MUA Spatial Grouping):\n', iRec, length(recordings));
    fprintf('Subject: %s, Session: %s, AcqID: %s, Sorter: %s\n', ...
        Subject, Session, AcquisitionsID, SorterName);
    fprintf('============================================================\n');
    
    % try
        %% 1. Load ephys and nidaq data info
        fprintf('Step 1: Loading experiment info...\n');
        [ExpInfo, ephysDir, nidqDir, ~] = getExpInfo(BaseDir, Subject, Session, AcquisitionsID);
        
        %% 2. Load units
        fprintf('Step 2: Loading units...\n');
        units = createClusterTable(ephysDir, SorterName);
        units = table2struct(units);
        
        %% 2.5 Group units into MUA by Shank and Depth
        fprintf('Step 2.5: Grouping units into MUA spatial bins (Bin size = %d um)...\n', mua_depth_bin_size);
        
        % Extract all depths and shanks, and FORCE them to be column vectors
        all_depths = arrayfun(@(x) x.location(2), units);
        all_depths = all_depths(:); 
        
        all_shanks = [units.shank];
        all_shanks = all_shanks(:);
        
        % Define depth edges based on the min/max depths of found units
        min_depth = floor(min(all_depths) / mua_depth_bin_size) * mua_depth_bin_size;
        max_depth = ceil(max(all_depths) / mua_depth_bin_size) * mua_depth_bin_size;
        depth_edges = min_depth : mua_depth_bin_size : max_depth;
        num_depth_bins = length(depth_edges) - 1;
        unique_shanks = unique(all_shanks)'; % Transpose for clean looping
        
        % Initialize the new MUA structure
        mua_sites = struct('shank', {}, 'depth_center', {}, 'spike_times', {}, 'num_units_pooled', {});
        site_idx = 1;
        
        for s = unique_shanks
            for d = 1:num_depth_bins
                d_start = depth_edges(d);
                d_end = depth_edges(d+1);
                d_center = (d_start + d_end) / 2;
                
                % Find all units that fall into this specific shank and depth bin
                valid_units = find(all_shanks == s & all_depths >= d_start & all_depths < d_end);
                
                % Force valid_units to be a row vector so the 'for' loop iterates correctly
                valid_units = valid_units(:)';
                
                % Only create a site if there is actually neural activity here
                if ~isempty(valid_units)
                    % Pool and sort spike times
                    pooled_spikes = [];
                    for u = valid_units
                        % Force column vector concatenation to prevent dimension mismatch
                        pooled_spikes = [pooled_spikes; units(u).spike_times(:)]; 
                    end
                    pooled_spikes = sort(pooled_spikes);
                    
                    mua_sites(site_idx).shank = s;
                    mua_sites(site_idx).depth_center = d_center;
                    mua_sites(site_idx).spike_times = pooled_spikes;
                    mua_sites(site_idx).num_units_pooled = length(valid_units);
                    site_idx = site_idx + 1;
                end
            end
        end
        fprintf('  -> Pooled %d single units into %d spatial MUA sites.\n', length(units), length(mua_sites));

        %% 3. Generate stimulus matrix from bonsai binary file
        fprintf('Step 3: Generating stimulus matrix...\n');
        binFileDir = fullfile(BaseDir, Subject, 'ephys', Session, ['Processed_', AcquisitionsID], 'bonsai_files');
        files = dir(fullfile(binFileDir, '*SparseNoise*.bin'));
        
        if isempty(files)
            error('No SparseNoise*.bin file found.');
        elseif length(files) > 1
            warning('Multiple matching files found. Using the first one: %s', files(1).name);
        end
        
        fileID = fopen(fullfile(files(1).folder, files(1).name));
        thisBinFile = fread(fileID);
        fclose(fileID);
        
        % Translate stimulus into -1:1 scale
        stim_matrix = zeros(1, length(thisBinFile));
        stim_matrix(thisBinFile == 0) = -1;
        stim_matrix(thisBinFile == 255) = 1;
        stim_matrix(thisBinFile == 128) = 0;
        
        % Reshape stimulus log into an NxM grid over time
        stim_matrix = reshape(stim_matrix, [grid_size(1), grid_size(2), length(thisBinFile) / (grid_size(1) * grid_size(2))]);
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
        pd_min = prctile(thisPD, 5); 
        pd_max = prctile(thisPD, 95); 
        pd_range = pd_max - pd_min;
        
        thresh_up = pd_min + (0.6 * pd_range);
        thresh_down = pd_min + (0.4 * pd_range);
        
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
        
       %% 5. Bin spike times relative to stimulus frames (High Temporal Resolution)
        fprintf('Step 5: Binning spikes with high temporal resolution...\n');
        
        % --- DECOUPLED NEURAL TIMING ---
        bin_size = 0.020; % 10 ms bins (High resolution for spike dynamics)
        max_response_time = 0.200; % Look up to 200ms post-stimulus
        AnalysisTimeWindow = [0, max_response_time]; 
        
        % This will create exactly 20 temporal bins per trial (0.200 / 0.010 = 20)
        num_temporal_bins = round(max_response_time / bin_size);
        
        resps = zeros(length(mua_sites), num_temporal_bins, length(stim_times));
        tic
        for isite = 1:length(mua_sites)
            spikes_this_site = mua_sites(isite).spike_times;
            % Sanity check to prevent passing empty spikes to PSTH
            if ~isempty(spikes_this_site)
                [~, ~, ~, ~, ~, binnedArray] = psthAndBA(spikes_this_site, stim_times, AnalysisTimeWindow, bin_size);
                resps(isite, :, :) = binnedArray';
            end
        end
        toc
        
        %% 6. Fit receptive fields
        fprintf('Step 6: Fitting receptive fields...\n');
        
        sn_options.framesToShow = 1:num_temporal_bins; 
        sn_options.plotflag = 0; 
        sn_options.grid_size = grid_size;
        
        rf_start_time = 0.040; 
        rf_end_time = 0.120;
        
        start_bin = round(rf_start_time / bin_size) + 1;
        end_bin = round(rf_end_time / bin_size);
        
        rf_names = {'linear', 'black', 'white', 'contrast'};
        
        tic
        parfor isite = 1:length(mua_sites)
            MUA_data = squeeze(resps(isite, :, :))'; 
            
            % Pass the original stim_matrix directly since we strictly validated its size in 4.5
            tempRF = sparseNoiseAnalysis(stim_matrix, MUA_data, [], [], sn_options);
            mua_sites(isite).allRF = tempRF; 
            
            for type_idx = 1:4
                this_rf_time = tempRF(:, :, :, type_idx);
                avg_rf = mean(this_rf_time(:, :, start_bin:end_bin), 3); 
                
                temp_up = imresize(avg_rf, scal_f);
                smth_rf = imgaussfilt(temp_up, sigma);
                
                norm_rf = smth_rf ./ max(max(abs(smth_rf(:))), 1e-9); 
                
                name = rf_names{type_idx};
                mua_sites(isite).([name, 'RF_raw']) = avg_rf;
                mua_sites(isite).([name, 'RF_smth']) = smth_rf;
                mua_sites(isite).([name, 'RF_norm']) = norm_rf;
            end
            
            mua_sites(isite).RF = mua_sites(isite).contrastRF_norm; 
            mua_sites(isite).smthRF = mua_sites(isite).contrastRF_smth;
        end
        toc
        
        %% 7. Save results
        fprintf('Step 7: Saving results...\n');
        fname = fullfile('D:\Code\ephysAnalysis2024\SparseNoiseAnalysis', ...
            [Subject, '_', Session, '_', AcquisitionsID, '_MUA_RFs.mat']);
        
        try
            save(fname, 'mua_sites');
        catch
            save(fname, 'mua_sites', '-v7.3');
        end
        fprintf('Successfully saved results to: %s\n', fname);
        
    % catch ME
    %     % --- Error Handling ---
    %     fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    %     fprintf('   ERROR processing Subject %s, Session %s, AcqID %s\n', ...
    %         Subject, Session, AcquisitionsID);
    %     fprintf('   Error message: %s\n', ME.message);
    %     fprintf('   Error occurred in file %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    %     fprintf('   Skipping to the next recording...\n');
    %     fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    %     continue; 
    % end
end
fprintf('\n\nBatch processing complete.\n');