%% STA filters from white noise stimulus

Subject = 'M24019';
Session = '20240701';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged','sorter_output');

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

clear NidqBin

%% get white noise stim

stim_idx = find(strcmp(ExpInfo.Stimulus, 'WhiteNoise'));
if stim_idx == height(ExpInfo)
    nidq_idx = ExpInfo.smp_nidq(stim_idx):numel(NidqTime);
else
    nidq_idx = ExpInfo.smp_nidq(stim_idx):(ExpInfo.smp_nidq(stim_idx+1)-1);
end

%% get instantaneous frequency of PWM sign
threshold = 600;

t = NidqTime(nidq_idx);
this_pd = PD(nidq_idx);

pwm = this_pd;

s=pwm>threshold;%pwm is your signal, threshold is when it is considered on or off
d=diff(s);
up_idx = find(d==1)+1; % +1 accouunts for diff
up_times = t(up_idx);
down_idx = find(d==-1)+1;
down_times = t(down_idx);

for iup = 1:numel(up_times)-1
    dc(iup) = (down_times(iup)-up_times(iup))/(up_times(iup+1)-up_times(iup));
end

% samples between up_idx

upd = diff(up_idx);

dct = repelem(dc,upd);
dct = [zeros(1, up_idx(1)), dct];
dct = [dct, zeros(1,numel(t)-numel(dct))];

down_times_end = down_times(end);

clear up_idx up_times down_idx down_times s pwm upd d dc PD NidqTime
%% load units

Subject = 'M24019';
Session = '20240701';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
SorterName = 'kilosort3_merged';
units = createClusterTable(Subject, Session, BaseDir, SorterName);
units = table2struct(units);

%%
timeWindow = 0.2;
nSamps = ceil(timeWindow * niSampRate);

startTime=tic;
for iunit = 1:numel(units)
iunit
spikeTimes = units(iunit).spike_times;
spikeTimes = spikeTimes(spikeTimes>t(1)+timeWindow);
spikeTimes = spikeTimes(spikeTimes<down_times_end);



STA = nan(numel(spikeTimes),nSamps);
runSpeedSpike = nan(1,numel(spikeTimes));

tic
for ispike = 1:numel(spikeTimes)
    spike_idx = find(t>spikeTimes(ispike),1,'first');
    STA(ispike,:) = dct(spike_idx-(nSamps-1):spike_idx);
    wheel_idx = find(wheel.time>spikeTimes(ispike),1,"first");
    runSpeedSpike(ispike) = wheel.smthSpeed(wheel_idx);
end
toc

sta_units2(iunit).STA = STA;
sta_units2(iunit).runSpeedSpike = runSpeedSpike;

end
stopTime = toc(startTime)

%% plot

statThresh = 1;
runThresh = 3;

for iunit = 1:10;%1:numel(sta_units2)

STA = sta_units2(iunit).STA;
runSpeedSpike = sta_units2(iunit).runSpeedSpike;

figure
title(['White noise STA for unit: ', num2str(iunit)])
xvals = (1:nSamps)./niSampRate;

shadedErrorBar(xvals, nanmean(STA(runSpeedSpike<statThresh,:),1),...
    nansem(STA(runSpeedSpike<statThresh,:),1))
shadedErrorBar(xvals, nanmean(STA(runSpeedSpike>runThresh,:),1),...
    nansem(STA(runSpeedSpike>runThresh,:),1),...
    'lineProps','r')
ax = gca;
ax.XTick = (0:0.05:0.2);
ax.XTickLabel = -200:50:0;
ylabel('PWM duty cycle')
xlabel('Time (ms)')
xlim([0 0.2])
defaultAxesProperties(gca, true)
pause
close
end