%% Basic analysis pipeline

Subject = 'M24065';
Session = '20250728';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort4','sorter_output');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session);
clusterTable = createClusterTable(kilosortDir);
units = table2struct(clusterTable);

%% get wheel data from nidq
samplingFreq=100; %hz
wheelChans = [3,4]; % NIDQ channels for wheel signal
wheel = getWheelPos(nidqDir,wheelChans);
wheel = processNidqWheel(wheel,samplingFreq,'gaussian',0.175); % window size in s

%% get photodiode data from t prime -> works for normal stim
photodiode.up_times = table2array(readtable(fullfile(nidqDir,'photodiode_up_tprime.txt')));
photodiode.down_times = table2array(readtable(fullfile(nidqDir,'photodiode_down_tprime.txt')));

%% generate trial structs

stim_index0 = [2,4];
stim_index50 = [3,5];

dmtrials_0 = createTrialsStruct(ExpInfo,stim_index0,photodiode,wheel);
stat_idx = [dmtrials_0.meanWheelSpeed]<=1;
run_idx = [dmtrials_0.meanWheelSpeed]>5;
[dmtrials_0.runFlag] = deal([nan]);
[dmtrials_0(stat_idx).runFlag] = deal([0]);
[dmtrials_0(run_idx).runFlag] = deal([1]);

dmtrials_50 = createTrialsStruct(ExpInfo,stim_index50,photodiode,wheel);
stat_idx = [dmtrials_50.meanWheelSpeed]<=1;
run_idx = [dmtrials_50.meanWheelSpeed]>5;
[dmtrials_50.runFlag] = deal([nan]);
[dmtrials_50(stat_idx).runFlag] = deal([0]);
[dmtrials_50(run_idx).runFlag] = deal([1]);

%% get spike counts
% binned spike count options
scopts.intervalStart = 0;
scopts.intervalEnd = 1.2;
scopts.binSpacing=1.2;
scopts.uniqueVals = {[0 16 32 64 128 256],[0, 1]};

[anUnits, cond] = getBinnedSpikeCounts(dmtrials_0, units, {'VelX1','runFlag'}, scopts);

for iunit = 1:numel(units)
    units(iunit).allSpikes_0 = anUnits(iunit).allSpikes;
end

[anUnits, cond] = getBinnedSpikeCounts(dmtrials_50, units, {'VelX1','runFlag'}, scopts);

for iunit = 1:numel(units)
    units(iunit).allSpikes_50 = anUnits(iunit).allSpikes;
end


%% downsample and get tuning strength

%% downsample to min # trials available
% set so equal # trials per condition

minTrial = min([min(min(cellfun(@(x) size(x,2), units(1).allSpikes_0))), min(min(cellfun(@(x) size(x,2), units(1).allSpikes_50)))]);

% units = units;

for iunit = 1:numel(units)
    units(iunit).allSpikes_0_ds = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes_0,'UniformOutput', false);
    units(iunit).allSpikes_50_ds = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes_50,'UniformOutput', false);

end

%% get cross-validated R2 (tuning strength) 

% options
r2opts.nPerms = 10;
r2opts.randFlag = 1;
r2opts.validMeans = 1;
r2opts.kval = 3;
r2opts.nShuffle = 100;

%get the tuning strength and significance
for iunit = 1:numel(units)

    % 0%
    %stationary trials
    sca = units(iunit).allSpikes_0_ds(:,1)';
    [units(iunit).statR2_0, units(iunit).statR2_pval_0] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

    %locomotion trials
    sca = units(iunit).allSpikes_0_ds(:,2)';
    [units(iunit).runR2_0, units(iunit).runR2_pval_0] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

     % 50%
    %stationary trials
    sca = units(iunit).allSpikes_50_ds(:,1)';
    [units(iunit).statR2_50, units(iunit).statR2_pval_50] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

    %locomotion trials
    sca = units(iunit).allSpikes_50_ds(:,2)';
    [units(iunit).runR2_50, units(iunit).runR2_pval_50] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

end

figure,
subplot(221)
hold on
plot([-1 1],[-1 1],'k')
plot([0.1 0.1], [-1 1],'k:')
plot([-1 1], [0.1 0.1],'k:')
plot([units.statR2_0],[units.statR2_50],'k.')
xlabel('0%'), ylabel('50%')
title('stat')

subplot(222)
hold on
plot([-1 1],[-1 1],'k')
plot([0.1 0.1], [-1 1],'k:')
plot([-1 1], [0.1 0.1],'k:')
plot([units.runR2_0],[units.runR2_50],'k.')
xlabel('0%'), ylabel('50%')
title('run')

subplot(223)
hold on
plot([-1 1],[-1 1],'k')
plot([0.1 0.1], [-1 1],'k:')
plot([-1 1], [0.1 0.1],'k:')
plot([units.statR2_0],[units.runR2_0],'k.')
xlabel('stat'), ylabel('run')
title('0%')

subplot(224)
hold on
plot([-1 1],[-1 1],'k')
plot([0.1 0.1], [-1 1],'k:')
plot([-1 1], [0.1 0.1],'k:')
plot([units.statR2_50],[units.runR2_50],'k.')
xlabel('stat'), ylabel('run')
title('50%')


%%
for iunit = 1:numel(units)
    ax(1) = subplot(121)
    shadedErrorBar(1:6, cellfun(@mean, units(iunit).allSpikes_0(:,1)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_0(:,1)));
    shadedErrorBar(1:6, cellfun(@mean, units(iunit).allSpikes_0(:,2)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_0(:,2)),'lineProps','r');
  ax(2) = subplot(122)
    shadedErrorBar(1:6, cellfun(@mean, units(iunit).allSpikes_50(:,1)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_50(:,1)));
    shadedErrorBar(1:6, cellfun(@mean, units(iunit).allSpikes_50(:,2)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_50(:,2)),'lineProps','r');
    linkaxes(ax)
    pause
    close
end