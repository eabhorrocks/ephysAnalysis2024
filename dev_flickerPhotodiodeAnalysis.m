%% flicker photodiode analysis

Subject = 'M24019';
Session = '20240718';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort4','sorter_output');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session);
% clusterTable = createClusterTable(kilosortDir);
% units = table2struct(clusterTable);

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

allStimTimes = [];
allStimTimesImec0=[];

for stim_index = [3,5]

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


%% get photodiode and wheel data
% quick hack for photodidoe
photodiode.up_times = allStimTimesImec0;
photodiode.down_times = photodiode.up_times+2;




%%
allFlickerTrials = createTrialsStruct(ExpInfo,[3,5],photodiode,wheel);

stat_idx = [allFlickerTrials.meanWheelSpeed]<=2;
run_idx = [allFlickerTrials.meanWheelSpeed]>5;
[allFlickerTrials.runFlag] = deal([nan]);
[allFlickerTrials(stat_idx).runFlag] = deal([0]);
[allFlickerTrials(run_idx).runFlag] = deal([1]);


%% get units
flickFreqs = unique([allFlickerTrials.StimulusFrequency]);

options.intervalEnd=2;
options.binSpacing = 2;
options.uniqueVals = {[flickFreqs], [0 1]};
[anUnits, cond] = getBinnedSpikeCounts(allFlickerTrials, units, {'StimulusFrequency', 'runFlag'}, options);


for iunit = 1:numel(units)
    units(iunit).allSpikes_0b = anUnits(iunit).allSpikes;
end

%%


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
% timeforspikes = [ExpInfo.sec_nidq(2), ExpInfo.sec_nidq(3)];
% 
% temp_spikes = cat(1,units.spike_times);
% 
% [binSpikes,edges] = histcounts(temp_spikes, timeforspikes(1):0.01:timeforspikes(2));
% smthSpikes = smoothdata(binSpikes,'gaussian',17.5);
% binTime = movmean(edges,2,'Endpoints','discard');
% 
% figure
% ax(1)=subplot(211);
% plot(thisTime,PD_binary)
% ax(2)=subplot(212);
% plot(binTime,smthSpikes)
% linkaxes(ax,'x');

