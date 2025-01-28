%% flicker photodiode analysis

Subject = 'M24077';
Session = '20241219';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session);
clusterTable = createClusterTable(kilosortDir);
units = table2struct(clusterTable);

% nidqDir = 'Z:\ibn-vision\DATA\SUBJECTS\M24019\ephys\20240718\nidq_processed';

%% get wheel data from nidq
samplingFreq=100; %hz
wheelChans = [3,4]; % NIDQ channels for wheel signal
wheel = getWheelPos(nidqDir,wheelChans);
wheel = processNidqWheel(wheel,samplingFreq,'gaussian',0.175); % window size in s

%%
nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

PD = NidqBin(2,:);

niSampRate = str2double(NidqMeta.niSampRate);
NidqTime = (1:size(NidqBin,2))./niSampRate;

%% get sync times for imec0 and nidq streams
nidq_syncTimes_imec0 = readmatrix(fullfile(nidqDir,'nidq_sync_tprime.txt'));
nidq_syncData = NidqBin(9,:);
syncIdx = find(diff(nidq_syncData)==8)+1;
nidq_syncTimes_nidq = NidqTime(syncIdx);


%% get stimulus start times for flicker

stim_idx_to_use = [1];

allStimTimes = [];
allStimTimesImec0=[];

for stim_index = stim_idx_to_use

    stimIdx = stim_index;

    if stimIdx<height(ExpInfo)
        idx = [ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx+1)];
    else
        idx = [ExpInfo.smp_nidq(stimIdx):numel(PD)];
    end

    % get PD signal and timestamps
    thisPD = PD(idx);
    thisTime = NidqTime(idx);
    stimTimeIdx = idx;

    % convert to a digital signal
    PD_binary = thisPD;
    threshold = 600;
    PD_binary(PD_binary<threshold)=0;
    PD_binary(PD_binary>=threshold)=1;

    % find trial start times
    idx = find(diff(PD_binary)==1); % indexes of where PD goes up
    idx = idx+1;
    diffTimes = thisTime(idx); % time stamps of where PD goes up
    diffDiffTimes = diff(diffTimes); % time between PD going up.
    idx2 = find(diffDiffTimes>1.5); % For new trials its ~2, i.e. >1.5s
    idx2=idx2+1;
    trialStartIdx = idx(idx2); % convert back to original PD index
    trialStartIdx = [idx(1), trialStartIdx]; % add the first trial
    stimTimes = thisTime(trialStartIdx); % get the PD timestamp for trial start idx

    timeFrame = [NidqTime(stimTimeIdx(1)), NidqTime(stimTimeIdx(end))-5];

    % convert stim times to imec0 timestamp using interp1
    % (maybe this should be converted to tprime)
    this_nidqsync_nidq = nidq_syncTimes_nidq(nidq_syncTimes_nidq >= timeFrame(1) & nidq_syncTimes_nidq <= timeFrame(2));
    this_nidqsync_imec0 = nidq_syncTimes_imec0(nidq_syncTimes_imec0 >= timeFrame(1) & nidq_syncTimes_imec0 <= timeFrame(2));
    this_nidqsync_nidq=this_nidqsync_nidq(:); this_nidqsync_imec0=this_nidqsync_imec0(:);


    while numel(this_nidqsync_imec0)>numel(this_nidqsync_nidq)
        dv = this_nidqsync_imec0(1:numel(this_nidqsync_nidq))-this_nidqsync_nidq;
        idx = find(abs(dv)>0.5,1,'first');
        this_nidqsync_imec0(idx)=[];
    end

    while numel(this_nidqsync_nidq)>numel(this_nidqsync_imec0)
        dv = this_nidqsync_nidq(1:numel(this_nidqsync_imec0))-this_nidqsync_imec0;
        idx = find(abs(dv)>0.5,1,'first');
        this_nidqsync_nidq(idx)=[];
    end


    stimTimes_imec0 = interp1(this_nidqsync_nidq, this_nidqsync_imec0, stimTimes);

    allStimTimes = cat(2,allStimTimes, stimTimes);
    allStimTimesImec0 = cat(2,allStimTimesImec0, stimTimes_imec0);


end

% plot binary PD signal with trial onset times
% figure, hold on
% plot(thisTime,PD_binary);
% plot(stimTimes, ones(size(stimTimes)), 'r*')
% xlim([1500 1520])
%


% get photodiode data
% quick hack for photodidoe
photodiode.up_times = allStimTimesImec0;
photodiode.down_times = photodiode.up_times+2;




%
allFlickerTrials = createTrialsStruct(ExpInfo,stim_idx_to_use,photodiode,wheel);

stat_idx = [allFlickerTrials.meanWheelSpeed]<=1;
run_idx = [allFlickerTrials.meanWheelSpeed]>3;
[allFlickerTrials.runFlag] = deal([nan]);
[allFlickerTrials(stat_idx).runFlag] = deal([0]);
[allFlickerTrials(run_idx).runFlag] = deal([1]);




%% LFP analysis

lfpDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0');
lfpMetaFile =  dir(fullfile(lfpDir,'*.meta'));
lfpBinFile =  dir(fullfile(lfpDir,'*.bin'));

lfpMeta = ReadMeta(fullfile(lfpMetaFile.folder, lfpMetaFile.name));
lfpSampRate = str2double(lfpMeta.imSampRate);
nSamples = str2double(lfpMeta.fileTimeSecs)*lfpSampRate;
lfpTime = (1:nSamples)/lfpSampRate;


%% get lfp for each trial

params.Fs=2500;
params.fpass = [1 150];
params.trialave = 0;
params.pad = 0;
params.tapers = [3 3];
params.err = [1 0.05];

chs = [377,327,189, 214];
for ichannel = 1%:4
    channel = chs(ichannel)
    preTime = 0;
    postTime = 0;

    for itrial = 1:numel(allFlickerTrials)
        start_idx = find(lfpTime<=allFlickerTrials(itrial).start_time-preTime,1,'last');
        end_idx = find(lfpTime>=allFlickerTrials(itrial).stop_time+postTime,1,'first');
        nSamp = end_idx-start_idx;

        dataArray = ReadBin(start_idx,nSamp, lfpMeta, lfpBinFile.name, lfpBinFile.folder);
        dataArrayCh = dataArray(channel,:);


        allFlickerTrials(itrial).LFP = dataArrayCh;


        [S,f,err] = mtspectrumc( dataArrayCh(:), params );

        allFlickerTrials(itrial).S = S;
        allFlickerTrials(itrial).normS = pow2db(S);%close pow2db(mean(S,2)/mean(mean(S,2)));

    end



    %

    figure, hold on
    uniqueFreqs = unique([allFlickerTrials.StimulusFrequency]);

    for ifreq = 1:numel(uniqueFreqs)
        subplot(2,5,ifreq), hold on
        idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
            [allFlickerTrials.runFlag]==0);

        allVals = cat(2,[allFlickerTrials(idx).normS])';
        plot(f,mean(allVals),'k')

        idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
            [allFlickerTrials.runFlag]==1);

        allVals = cat(2,[allFlickerTrials(idx).normS])';
        plot(f,mean(allVals),'r')

        title(uniqueFreqs(ifreq))
        xlabel('Frequency (Hz)')
        ylabel('normalised power')
    end
end
% legend(uniqueFreqs)


%% subtract 90hz as baseline

ifreq=10;
idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
    [allFlickerTrials.runFlag]==0);

allVals = cat(2,[allFlickerTrials(idx).normS])';
statBaseline = mean(allVals);

idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
    [allFlickerTrials.runFlag]==1);

allVals = cat(2,[allFlickerTrials(idx).normS])';
runBaseline = mean(allVals);

for ifreq = 1:numel(uniqueFreqs)
        subplot(2,5,ifreq), hold on
        idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
            [allFlickerTrials.runFlag]==0);

        allVals = cat(2,[allFlickerTrials(idx).normS])';
        plot(f,mean(allVals)-statBaseline,'k')

        idx = find([allFlickerTrials.StimulusFrequency]==uniqueFreqs(ifreq) &...
            [allFlickerTrials.runFlag]==1);

        allVals = cat(2,[allFlickerTrials(idx).normS])';
        plot(f,mean(allVals)-runBaseline,'r')

        title(uniqueFreqs(ifreq))
        xlabel('Frequency (Hz)')
        ylabel('normalised power')
        ylim([-4 13])
    end


%% get units


% Subject = 'M24019';
% Session = '20240718';
% BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
SorterName = 'kilosort3_merged';
units = createClusterTable(Subject, Session, BaseDir, SorterName);
units = table2struct(units);

%% get spike counts

flickFreqs = unique([allFlickerTrials.StimulusFrequency]);

options.intervalEnd=2;
options.binSpacing = 2;
options.uniqueVals = {[flickFreqs], [0 1]};
[anUnits, cond] = getBinnedSpikeCounts(allFlickerTrials, units, {'StimulusFrequency', 'runFlag'}, options);


for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

%% plot tuning curves, one unit at a time


figure
for iunit = 1:numel(units)
    ax(1) = subplot(121);

    % get custom ylims
    % minVal = floor(min([min(cellfun(@(x) mean(x,2), units(iunit).allSpikes_0b)-cellfun(@(x) sem(x,2), units(iunit).allSpikes_0b),[],'all'), ...
    %     min(cellfun(@(x) mean(x,2), units(iunit).allSpikes_50b)-cellfun(@(x) sem(x,2), units(iunit).allSpikes_50b),[],'all')]));
    % maxVal = ceil(max([min(cellfun(@(x) mean(x,2), units(iunit).allSpikes_0b)+cellfun(@(x) sem(x,2), units(iunit).allSpikes_0b),[],'all'), ...
    %     max(cellfun(@(x) mean(x,2), units(iunit).allSpikes_50b)+cellfun(@(x) sem(x,2), units(iunit).allSpikes_50b),[],'all')]));


    hold off
    errorbar(flickFreqs, cellfun(@(x) mean(x,2), units(iunit).allSpikes_0b(:,1)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_0b(:,1)),'k')
    hold on
    errorbar(flickFreqs, cellfun(@(x) mean(x,2), units(iunit).allSpikes_0b(:,2)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_0b(:,2)),'r')
    title([num2str(iunit), ' low brightness'])
    xlabel('Flicker frequency')
    ylabel('Spike counts (2s duration)')
    defaultAxesProperties(gca, true)
    ax(2) = subplot(122);
    hold off
    errorbar(flickFreqs, cellfun(@(x) mean(x,2), units(iunit).allSpikes_50b(:,1)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_50b(:,1)),'k')
    hold on
    errorbar(flickFreqs, cellfun(@(x) mean(x,2), units(iunit).allSpikes_50b(:,2)), cellfun(@(x) sem(x,2), units(iunit).allSpikes_50b(:,2)),'r')
    % ylim([minVal, maxVal])
    title([num2str(iunit), ' high brightness'])
    xlabel('Flicker frequency')
    ylabel('Spike counts (2s duration)')
    defaultAxesProperties(gca, true)
    % linkaxes(ax,'y')
    pause
end



%%
%% get trial intervals for PSTH function

statTrials = allFlickerTrials([allFlickerTrials.runFlag]==0);
runTrials = allFlickerTrials([allFlickerTrials.runFlag]==1);

for ifreq = 1:numel(flickFreqs)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.StimulusFrequency]==flickFreqs(ifreq));
    stat_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)].*1000;

    stimTrials = runTrials([runTrials.StimulusFrequency]==flickFreqs(ifreq));
    run_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)].*1000;
end

%% generate PSTHs

% options for PSTH function
options.binWidth = 10;
options.preTime = 1000;
options.postTime = 1000;
options.smoothType = 'gaussian';
options.smoothWidth = 75;

options.getReliability = false;
options.kfold = 3; % set to number>nTrials to do Leave-one-out.
options.distType = 'correlation';
options.nReliPerms = 25;
options.nShuffle = 40;

options.plot = true;
options.cols = tab10; %options.cols(1,:) = [];
options.plotFRrange = [];
options.plot95CI = false;
figure
freq_idx = [2]
for iunit = 60%1:numel(units)
    iunit
    %[psth, blank, bmp] = makePSTH(spikeTimes,intervalsCell,blankIntervals,options);
    subplot(121)
    [units(iunit).stat_psth, ~, ~]            = makePSTH(units(iunit).spike_times*1000,stat_intervalsCell(freq_idx),[],options);
    subplot(122)
    [units(iunit).run_psth, ~, bmp]            = makePSTH(units(iunit).spike_times*1000,run_intervalsCell(freq_idx),[],options);
    pause
    close
end


%% get frequency spectrum of spiking responses

cols=inferno(5);
iunit = 4;

figure,
subplot(121),
hold on
for ifreq = 2
    plot(bmp, units(iunit).stat_psth(ifreq).psth,'Color',cols(ifreq,:))
end
subplot(122),
hold on
for ifreq = 2
    plot(bmp, units(iunit).run_psth(ifreq).psth,'Color',cols(ifreq,:))
end



