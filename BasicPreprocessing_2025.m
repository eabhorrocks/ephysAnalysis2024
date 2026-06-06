%% Basic Preprocessing
BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';
Subject = 'M25066';
Session = '20250728';
AcquisitionsID = '0';
SorterName = 'kilosort4';

%% Set up logging
% Create a timestamp string (e.g., '20250728_143005')
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Define the log filename
logFileName = sprintf('Log_%s_%s_%s.txt', Subject, Session, timestamp);

% Define where to save it (e.g., in the BaseDir or a specific 'Logs' folder)
% Here we try to put it in the session folder, but fallback to current folder if it doesn't exist yet
savePath = fullfile(BaseDir,Subject,'ephys',Session,['Processed_',AcquisitionsID],'preprocessed');

if ~isfolder(savePath)
    mkdir(savePath)
end

logFullPath = fullfile(savePath, logFileName);


% Start Logging
if exist(logFullPath, 'file'), delete(logFullPath); end % Clear if exists (unlikely with timestamp)
diary(logFullPath);
fprintf('--- LOG STARTED: %s ---\n', datestr(now));
fprintf('Saving log to: %s\n\n', logFullPath);

% try 
    %% load ephys and nidaq data info
    [ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);
    
    % Get 1hz sync pulse files
    nidq_1hz_file = dir(fullfile(nidqDir, '*nidq.xd_8_3_500.txt'));
    nidq_tprime_file = dir(fullfile(nidqDir, 'nidq_sync_tprime.txt'));
    
    % Build the full file paths
    nidq_1hz_file = fullfile(nidq_1hz_file.folder, nidq_1hz_file.name);
    nidq_tprime_file = fullfile(nidq_tprime_file.folder, nidq_tprime_file.name);
    
    % Specify NIDQ channels
    nidq_channels.asyncPulse = 1;
    nidq_channels.photodiode = 2;
    
    %% load units
    units = createClusterTable(ephysDir, SorterName);
    units = table2struct(units);
    
    %% Define all stimulus blocks to process and their associated properties
    stimBlocksToProcess(1).igkey = 0;
    stimBlocksToProcess(1).properties = {'speed'};
    stimBlocksToProcess(2).igkey = 1;
    stimBlocksToProcess(2).properties = {'coherence'};
    stimBlocksToProcess(3).igkey = 2;
    stimBlocksToProcess(3).properties = {'speed'};
    stimBlocksToProcess(4).igkey = 3;
    stimBlocksToProcess(4).properties = {'coherence'};
    
    % Initialize a cell array to hold the trial structs from each block
    all_trials_cell = cell(length(stimBlocksToProcess), 1);
    all_wheel_cell = cell(length(stimBlocksToProcess), 1);
    doPlot = false; % Disable plotting within the loop
    
    %% Loop over each defined stimulus block
    for i = 1:length(stimBlocksToProcess)
        
        current_igkey = stimBlocksToProcess(i).igkey;
        current_properties = stimBlocksToProcess(i).properties;
        
        fprintf('--- Processing block igkey: %d ---\n', current_igkey);
        
        % Get NIDQ data for this gX recording
        % Get minimally processed NIDQ data
        nidq = extractNidqData(ExpInfo, current_igkey, str2double(AcquisitionsID), nidqDir, nidq_channels, nidq_1hz_file, nidq_tprime_file, doPlot);
        
        % Generate trials struct from NIDQ and trialParams .csv file
        g_idx = find(ExpInfo.gKey == current_igkey);
        this_trialParamsCSVfile = ExpInfo.TrialParamsFiles{g_idx};
        trials = crateTrialsStruct_bonsai(nidq, this_trialParamsCSVfile);
        [trials.gkey]=deal(current_igkey);
        
        % Process the wheel data file
        wheel_tbl = readtable(ExpInfo.WheelFiles{g_idx});
        wheel_tbl.gkey = ones(height(wheel_tbl), 1) * current_igkey;
        wheel_tbl.dist = wheel2unit(wheel_tbl.Wheel,4096, 19.5);
        wheel_tbl.ddist_dt = [nan; diff(wheel_tbl.dist)./diff(wheel_tbl.ArduinoTime/1000)];
        wheel_tbl.ArduinoTime_mid = movmean(wheel_tbl.ArduinoTime/1000,2);
        wheel_tbl(1,:)=[];
        syncTimes_bonsai = unique(wheel_tbl.LastSyncPulseTime)/1000;
        wheel_tbl.imecTime = mapTimestampsUsingAsyncPulse(syncTimes_bonsai,nidq.asyncPulseTimes_imec,wheel_tbl.ArduinoTime_mid);
        
        % Smooth the wheel speed
        sigma = 0.035; % standard deviation of gaussian in seconds 
        wheel_tbl.smthSpeed = smoothDataWithGaussian(wheel_tbl.ddist_dt, wheel_tbl.imecTime, sigma);
        
        % Align wheel data with each trial
        preStimTime = 0.5;
        postStimTime = 0.5;
        trials = alignDataToTrials(trials, wheel_tbl.imecTime, preStimTime, postStimTime,...
            'wheelSpeed', wheel_tbl.smthSpeed);
        
        % Inject the block-specific properties into every trial
        % 'deal' assigns the same value to a field in every element of the struct array
        if ~isempty(trials)
            [trials.properties] = deal(current_properties);
        end
        
        % Store the completed trials struct for this block
        all_trials_cell{i} = trials;
        all_wheel_cell{i} = wheel_tbl;
        
        fprintf('--- Finished block %d. Found %d matched trials. ---\n\n', current_igkey, length(trials));
        clear nidq trials wheel_tbl g_idx
    end

    %% save data
    saveName = fullfile(savePath, [Subject,'_',Session,'.mat']);
    % save()
    fprintf('Processing Complete.\n');

% catch ME
%     % If an error occurs, capture it in the log before closing
%     fprintf('\n!!! ERROR OCCURRED !!!\n');
%     fprintf('%s\n', ME.message);
%     fprintf('Error in: %s (Line %d)\n', ME.stack(1).name, ME.stack(1).line);
% end


%% --- STOP LOGGING ---
fprintf('--- LOG ENDED: %s ---\n', datestr(now));
diary off;
