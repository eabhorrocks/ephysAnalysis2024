function trials = createTrialsStruct(ExpInfo,stimIdx,photodiode,wheel)

 trials = table;

for istim = 1:numel(stimIdx)

    thisStim=stimIdx(istim);

    pd_up_times = photodiode.up_times; %table2array(readtable(fullfile(nidqDir,'photodiode_up_tprime.txt')));
    pd_down_times = photodiode.down_times; %table2array(readtable(fullfile(nidqDir,'photodiode_down_tprime.txt')));
    startTime = ExpInfo.('sec_imap0')(thisStim);
    endTime = ExpInfo.('sec_imap0')(thisStim+1); % need to check for last file!!, toherwise use PD_down_times + X?

    temp_pd_up = pd_up_times(pd_up_times>=startTime & pd_up_times<endTime);
    temp_pd_down =  pd_down_times(pd_up_times>=startTime & pd_up_times<endTime);

    temp_trials = readtable(ExpInfo.TrialParamsFiles{thisStim});
    temp_trials.("start_time") = temp_pd_up(:);
    temp_trials.("stop_time") = temp_pd_down(:);

    % import wheel file for this stimulus index
    % thisWheelFile = ExpInfo.WheelFiles(thisStim);
    % smoothMethod='gaussian'; windowSize=175; %ms
    % wheel = processWheel2024Temp(thisWheelFile{:},smoothMethod,windowSize);

    wheel_trace = [];
    % offset = 
    for itrial = 1:height(temp_trials)
        [~,wheel_start_idx] = min(abs(wheel.time-(temp_trials.start_time(itrial)-0.2)));
        [~,wheel_end_idx] = min(abs(wheel.time-(temp_trials.stop_time(itrial)+0.2)));
        wheel_trace{itrial}(1,:) = wheel.smthSpeed(wheel_start_idx:wheel_end_idx);
        wheel_trace{itrial}(2,:) = wheel.time(wheel_start_idx:wheel_end_idx)-temp_trials.start_time(itrial);
    end

    meanWheelSpeeds = cellfun(@(x) mean(x(1,:),2), wheel_trace,'UniformOutput',false);
    temp_trials.WheelSpeedTrace= wheel_trace';
    temp_trials.meanWheelSpeed = meanWheelSpeeds(:);

    trials = cat(1,trials,temp_trials);
    clear temp_trials wheel_trace
end

trials = table2struct(trials);