% Batch processing script for receptive field mapping

% =========================================================================
% ==                  DEFINE YOUR RECORDINGS HERE                        ==
% =========================================================================
% Add each unique recording as an element in the 'recordings' struct array.
% You can add as many as you need.

recordings = struct([]); % Initialize the struct array

recordings(1).Subject = 'M26017';
recordings(1).Session = '20260604';
recordings(1).AcquisitionsID = '0';
recordings(1).SorterName = 'kilosort4';

for irecording = 2:8
    recordings(irecording) = recordings(1);
    recordings(irecording).AcquisitionsID = [num2str(irecording-1)];
end


% recordings(3).Subject = 'M25106';
% recordings(3).Session = '20250904';
% recordings(3).AcquisitionsID = '0';
% recordings(3).SorterName = 'kilosort4';
% =========================================================================


% Main directory for all subject data
BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';

% RF mapping parameters
grid_size = [8, 16]; % 15 deg, normal cortex
% grid_size = [12, 24]; % 10 deg
% grid_size = [16, 32]; % 7.5deg, normal lgn 

scal_f = 10; % Scaling factor for RF upsampling
sigma = 3;   % Sigma for Gaussian smoothing of the RF

% =========================================================================
% ==                      START BATCH PROCESSING LOOP                    ==
% =========================================================================

for iRec = 1:length(recordings)
    
    % --- Get parameters for the current recording ---
    Subject = recordings(iRec).Subject;
    Session = recordings(iRec).Session;
    AcquisitionsID = recordings(iRec).AcquisitionsID;
    SorterName = recordings(iRec).SorterName;
    
    fprintf('\n\n============================================================\n');
    fprintf('Processing recording %d of %d:\n', iRec, length(recordings));
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
        
        %% 4. Get photodiode signal for frame times
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
        
        % Threshold PD signal to find frame flips
        flipTime = 0.1; 
        threshold = 2000;
        PD_binary = thisPD > threshold;
        
        % Find upswings (black-to-white or white-to-black transitions)
        up_idx = find(diff(PD_binary) == 1) + 1;
        diffUpTimes = diff(thisTime(up_idx));
        upSwings = up_idx([true, diffUpTimes > flipTime]);

        % Find downswings
        down_idx = find(diff(PD_binary) == -1) + 1;
        diffDownTimes = diff(thisTime(down_idx));
        downSwings = down_idx([true, diffDownTimes > flipTime]);
        
        stim_times = sort(thisTime([upSwings, downSwings]));
        
        %% 5. Bin spike times relative to stimulus frames
        fprintf('Step 5: Binning spikes for %d units...\n', length(units));
        bin_size = 1/60; % 60 Hz stimulus rate
        AnalysisTimeWindow = [0, bin_size * 7]; % Analyze 7 frames post-stimulus
        
        resps = zeros(length(units), length(AnalysisTimeWindow(1):bin_size:AnalysisTimeWindow(2))-1, length(stim_times));

        tic
        for iunit = 1:length(units)
            spikes_this_cell = units(iunit).spike_times;
            [~, ~, ~, ~, ~, binnedArray] = psthAndBA(spikes_this_cell, stim_times, AnalysisTimeWindow, bin_size);
            resps(iunit, :, :) = binnedArray';
        end
        toc
        
        %% 6. Fit receptive fields
        fprintf('Step 6: Fitting receptive fields...\n');
        sn_options.framesToShow = 1:7; % Frames to include in analysis
        sn_options.plotflag = 0; % Set to 1 to plot RFs during computation
        sn_options.grid_size = grid_size;

        
        
        tic
        for iunit = 1:length(units)
            SUA_data = squeeze(resps(iunit, :, :))'; 
            SUA_data = SUA_data(1:size(stim_matrix, 3), :); % Ensure dimensions match
            
            tempRF = sparseNoiseAnalysis(stim_matrix, SUA_data, [], [], sn_options);
            units(iunit).allRF = tempRF;
            
            linearRFs = tempRF(:,:,:,4); % 1=linear, 2=black, 3=white, 4=contrast
            units(iunit).RF = mean(linearRFs(:,:,3:7), 3); % Avg frames 3-7
            
            temp = imresize(units(iunit).RF, scal_f);
            units(iunit).smthRF = imgaussfilt(temp, sigma);
            units(iunit).normRF = units(iunit).smthRF ./ abs(max(units(iunit).smthRF(:)));
        end
        toc
        
        %% 7. Save results
        fprintf('Step 7: Saving results...\n');
        fname = fullfile('D:\Code\ephysAnalysis2024\SparseNoiseAnalysis', ...
            [Subject, '_', Session, '_', AcquisitionsID, '_RFs.mat']);
        
        try
            save(fname, 'units');
        catch
            save(fname, 'units', '-v7.3');
        end
        fprintf('Successfully saved results to: %s\n', fname);

    % catch ME
    %     % --- Error Handling ---
    %     % If an error occurs, print it to the console and continue to the next recording
    %     fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    %     fprintf('   ERROR processing Subject %s, Session %s, AcqID %s\n', ...
    %         Subject, Session, AcquisitionsID);
    %     fprintf('   Error message: %s\n', ME.message);
    %     fprintf('   Error occurred in file %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    %     fprintf('   Skipping to the next recording...\n');
    %     fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    %     continue; % This command skips to the next iteration of the for-loop
    % end
end

fprintf('\n\nBatch processing complete.\n');
