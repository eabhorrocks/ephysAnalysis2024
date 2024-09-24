%% temp analysis of tuning curves


Subject = 'M24019';
Session = '20240719';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort4','sorter_output');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session);
clusterTable = createClusterTable(kilosortDir);
units = table2struct(clusterTable);


%% generate trial structs
 
%   0% trials

stimIdx = [2 4];
trials = table;

for istim = 1:numel(stimIdx)

    thisStim=stimIdx(istim);

    pd_up_times = table2array(readtable(fullfile(nidqDir,'photodiode_up_tprime.txt')));
    pd_down_times = table2array(readtable(fullfile(nidqDir,'photodiode_down_tprime.txt')));
    startTime = ExpInfo.('sec_imap0')(thisStim);
    endTime = ExpInfo.('sec_imap0')(thisStim+1);

    temp_pd_up = pd_up_times(pd_up_times>=startTime & pd_up_times<endTime);
    temp_pd_down =  pd_down_times(pd_up_times>=startTime & pd_up_times<endTime);

    temp_trials = readtable(ExpInfo.TrialParamsFiles{thisStim});
    temp_trials.("start_time") = temp_pd_up;
    temp_trials.("stop_time") = temp_pd_down;

    % import wheel file for this stimulus index
    thisWheelFile = ExpInfo.WheelFiles(thisStim);
    smoothMethod='gaussian'; windowSize=175; %ms
    wheel = processWheel2024Temp(thisWheelFile{:},smoothMethod,windowSize);

    wheel_trace = [];
    for itrial = 1:height(temp_trials)
        trialDuration = (temp_trials.stop_time(itrial)-temp_trials.start_time(itrial))*1000; % convert to ms
        wheel_start_idx = find(wheel.time>temp_trials.ArduinoTime(itrial),1,'first');
        wheel_end_idx = find(wheel.time<temp_trials.ArduinoTime(itrial) + trialDuration,1,'last');
        wheel_trace{itrial} = wheel.smthSpeed(wheel_start_idx:wheel_end_idx);
    end

    meanWheelSpeeds = cellfun(@mean, wheel_trace);
    temp_trials.wheel = wheel_trace(:);
    temp_trials.meanWheel = meanWheelSpeeds(:);

    trials = cat(1,trials,temp_trials);
    clear temp_trials wheel_trace
end

dmtrials_0 = table2struct(trials);
clear trials


%%
trials = cat(1,temp_trials);

dmtrials_0 = table2struct(trials);

%   50% trials
stimIdx = [3 5];

for istim = 1:numel(stimIdx)

    thisStim=stimIdx(istim);
    pd_up_times = table2array(readtable(fullfile(nidqDir,'photodiode_up_tprime.txt')));
    pd_down_times = table2array(readtable(fullfile(nidqDir,'photodiode_down_tprime.txt')));
    startTime = ExpInfo.('sec_imap0')(thisStim);
    endTime = ExpInfo.('sec_imap0')(thisStim+1);

    temp_pd_up = pd_up_times(pd_up_times>=startTime & pd_up_times<endTime);
    temp_pd_down =  pd_down_times(pd_up_times>=startTime & pd_up_times<endTime);


    temp_trials = readtable(ExpInfo.TrialParamsFiles{thisStim});
    temp_trials.("start_time") = temp_pd_up;
    temp_trials.("stop_time") = temp_pd_down;

end

trials = cat(1,temp_trials);

dmtrials_50 = table2struct(trials);

%%

% binned spike count options
scopts.intervalStart = 0;
scopts.intervalEnd = 1.2;
scopts.binSpacing=1.2;
scopts.uniqueVals = {[0 16 32 64 128 256]};

[anUnits, cond] = getBinnedSpikeCounts(dmtrials_0, units, {'VelX1'}, scopts);

for iunit = 1:numel(units)
    units(iunit).allSpikes_0 = anUnits(iunit).allSpikes;
    units(iunit).meantuning_0 = cellfun(@mean, anUnits(iunit).allSpikes);
    units(iunit).semtuning_0 = cellfun(@(x) sem(x,2), anUnits(iunit).allSpikes);
    % shadedErrorBar(1:6, units(iunit).meantuning, units(iunit).semtuning);
    % pause
    % close
end


%% 50% trials
stimIdx = [3 5];

for istim = 1:numel(stimIdx)

    thisStim=stimIdx(istim);
    pd_up_times = table2array(readtable(fullfile(nidqDir,'photodiode_up_tprime.txt')));
    pd_down_times = table2array(readtable(fullfile(nidqDir,'photodiode_down_tprime.txt')));
    startTime = ExpInfo.('sec_imap0')(thisStim);
    endTime = ExpInfo.('sec_imap0')(thisStim+1);

    temp_pd_up = pd_up_times(pd_up_times>=startTime & pd_up_times<endTime);
    temp_pd_down =  pd_down_times(pd_up_times>=startTime & pd_up_times<endTime);


    temp_trials = readtable(ExpInfo.TrialParamsFiles{thisStim});
    temp_trials.("start_time") = temp_pd_up;
    temp_trials.("stop_time") = temp_pd_down;

end

trials = cat(1,temp_trials);

dmtrials_50 = table2struct(trials);

[anUnits, cond] = getBinnedSpikeCounts(dmtrials, units, {'VelX1'}, scopts);

for iunit = 1:numel(units)
    units(iunit).allSpikes_50 = anUnits(iunit).allSpikes;
    units(iunit).meantuning_50 = cellfun(@mean, anUnits(iunit).allSpikes);
    units(iunit).semtuning_50 = cellfun(@(x) sem(x,2), anUnits(iunit).allSpikes);
    % shadedErrorBar(1:6, units(iunit).meantuning, units(iunit).semtuning);
    % pause
    % close
end


%%
for iunit = 1:numel(units)
    shadedErrorBar(1:6, units(iunit).meantuning_0, units(iunit).semtuning_0);
    shadedErrorBar(1:6, units(iunit).meantuning_50, units(iunit).semtuning_50,'lineProps','r');
    pause
    close
end