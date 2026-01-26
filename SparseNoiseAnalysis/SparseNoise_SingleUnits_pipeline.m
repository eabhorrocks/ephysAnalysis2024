%% Sparse noise analysis pipeline

BaseDir = 'Y:\ibn-vision\DATA\SUBJECTS';
Subject = 'M26001';
Session = '20260122';
AcquisitionsID = '5';
SorterName = 'kilosort4';

 grid_size = [8, 16]; % 15 deg, normal cortex
% grid_size = [12, 24]; % 10 deg
% grid_size = [16, 32]; % 7.5deg, normal lgn

%% load ephys and nidaq data info


[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);

%% load units
units = createClusterTable(ephysDir, SorterName);
units = table2struct(units);


%% generate stimulus matrix

% find sparse noise bin file and minimally process
binFileDir = fullfile(BaseDir,Subject,'ephys',Session,['Processed_',AcquisitionsID], 'bonsai_files');

files = dir(fullfile(binFileDir, '*SparseNoise*.bin'));

% Display found files and handle different cases
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

fileID=fopen(fullfile(files(1).folder,files(1).name)); % load the bin file

thisBinFile=fread(fileID);
fclose(fileID);

% Translate stimulus into -1:1 scale
stim_matrix = zeros(1,length(thisBinFile));
stim_matrix(thisBinFile==0)=-1;
stim_matrix(thisBinFile==255)=1;
stim_matrix(thisBinFile==128)=0;
% Make a NxM grid from the stimulus log
stim_matrix = reshape(stim_matrix, [grid_size(1),grid_size(2), length(thisBinFile)/grid_size(1)/grid_size(2)]);
stim_matrix = stim_matrix(:,:,1:end-1); % The last upswing should be ignored


%% get photodiode signal for frame times
nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

PD = NidqBin(2,:);

niSampRate = str2double(NidqMeta.niSampRate);
NidqTime = (1:size(NidqBin,2))./niSampRate;

% take photodiode signal relevant to sparse noise (exclude other stimuli)
stimIdx = find(strcmp(ExpInfo.Stimulus,'SparseNoise'))
if stimIdx<height(ExpInfo)
    idx = [ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx+1)];
else
    idx = [ExpInfo.smp_nidq(stimIdx):numel(NidqTime)];
end

thisPD = PD(idx);
thisTime = NidqTime(idx);

flipTime =0.1; % NB: add BUT NOT MORE THAN flipTime*2.5 or something

% threshold PD signal and get up and down times
% %convert to a digital signal
PD_binary = thisPD;
threshold = 2000;
PD_binary(PD_binary<threshold)=0;
PD_binary(PD_binary>=threshold)=1;

% find up swings
idx = find(diff(PD_binary)==1); % indexes of where PD goes up
idx = idx+1;
diffTimes = thisTime(idx); % time stamps of where PD goes up
diffDiffTimes = diff(diffTimes); % time between PD going up.
idx2 = find(diffDiffTimes>flipTime); % 
idx2=idx2+1;
upSwings = idx(idx2); % convert back to original PD index
upSwings = [idx(1), upSwings]; % add the first trial

% find down swings
idx = find(diff(PD_binary)==-1); % indexes of where PD goes up
idx = idx+1;
diffTimes = thisTime(idx); % time stamps of where PD goes up
diffDiffTimes = diff(diffTimes); % time between PD going up.
idx2 = find(diffDiffTimes>flipTime); %
idx2=idx2+1;
downSwings = idx(idx2); % convert back to original PD index
downSwings = [idx(1), downSwings];

stim_times = sort(thisTime([upSwings, downSwings]));

%% get binned spikes (should update to use own function)

sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear', 'black','white','contrast'};
sn_options.mapMethod = 'mean'; % fitlm mean
sn_options.framesToShow = 1:7;  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 0; % 1 if you want to see details of the receptive fields

bin_size = 1/60;
AnalysisTimeWindow = [0 1/60*7];
tic
for iunit = 1:length(units)
    spikes_this_cell = units(iunit).spike_times;

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, stim_times, AnalysisTimeWindow, bin_size);
    resps(iunit,:,:) = binnedArray';
end
toc

%% fit RFs

% check how much data needed
% stim_matrix = stim_matrix(:,:,3000);
% resps = resps(:,:,3000);

% resps = resps(:,:,1:3000);
% stim_matrix = stim_matrix(:,:,1:3000);

scal_f=10;
sigma=3;
tic
for iunit = 1:numel(units)
    % iunit
    SUA_data = squeeze(resps(iunit,:,:))'; % returns 1 x N time bins x nFrames;
    SUA_data = SUA_data(1:size(stim_matrix,3),:); % sometimes theres an extra element
    tempRF = sparseNoiseAnalysis(stim_matrix,SUA_data,[],[],sn_options);
    units(iunit).allRF = tempRF;

    linearRFs =tempRF(:,:,:,4); % 1 = linear RFs, 4 = contrast RFs
    units(iunit).RF = mean(linearRFs(:,:,3:7),3);
    temp = imresize(units(iunit).RF,scal_f);
    units(iunit).smthRF = imgaussfilt(temp,sigma);
    units(iunit).normRF = units(iunit).smthRF./abs(max(units(iunit).smthRF(:)));

end
toc


%% save

fname = [Subject, '_', Session, '_', AcquisitionsID, '_RFs.mat'];
try
    save(fname, 'units');
catch
    save(fname, 'units', '-v7.3')
end


%% filter units if desired (e.g. for diff areas)

allocs = cat(1,units.location);
filtUnits=units;
% filtUnits = units(allocs(:,2)>4000);

%% plot average RFs by shank

figure
for ishank = 0:3
    subplot(2,2,ishank+1)
    idx = find([filtUnits.shank]==ishank);
    rf = nanmean(cat(3,filtUnits(idx).normRF),3);
    if ~isempty(rf)
        imagesc(flipud(rf))
        xlabel('Azimuth')
        ylabel('Elevation')
        xticks(linspace(1,size(rf,2),13))
        xticklabels(-120:20:120)
        yticks(linspace(1,size(rf,1),7))
        yticklabels(-60:20:60)
        title(sprintf('Shank %i Receptive field',ishank))
        % colorbar
        % colormap(1-gray)
        axis xy
    else
        warning('No units found on shank %d', ishank);
    end
end

%% plot RF by shank and depth
%split depth into
figure

nDepthSplits=10;
maxVal = max(allocs(:,2));
minVal = min(allocs(:,2));
edges = linspace(minVal,maxVal,nDepthSplits+1);
depth_disc = discretize(allocs(:,2),edges)';

iplot=0;
for idepth = nDepthSplits:-1:1

    for ishank=1:4

        iplot=iplot+1;
        try
            subplot(nDepthSplits,4, iplot)
            idx = find([units.shank]==ishank & depth_disc==idepth);
            rf = nanmean(cat(3,units(idx).normRF),3);
            imagesc(flipud(rf))
            xlabel('Azimuth')
            ylabel('Elevation')
            xticks(linspace(1,size(rf,2),13))
            xticklabels(-120:20:120)
            yticks(linspace(1,size(rf,1),7))
            yticklabels(-60:20:60)
            title(sprintf('Shank %i, Depth %i',ishank,idepth))
            % colorbar
            % colormap(1-gray)
            axis xy
        catch
        end
    end
end




%% plot one at a time

idx = find([filtUnits.firing_rate]>3);

% for iunit =186% 1:numel(idx)
%
%     thisMap_s = filtUnits(idx(iunit)).smthRF;
%     imagesc(flipud(thisMap_s))
%     xlabel('Azimuth')
%     ylabel('Elevation')
%     xticks(linspace(1,size(thisMap_s,2),13))
%     xticklabels(-120:20:120)
%     yticks(linspace(1,size(thisMap_s,1),7))
%     yticklabels(-60:20:60)
%     title(sprintf('Unit: %i, Shank: %i', idx(iunit), filtUnits(idx(iunit)).shank))
%     colorbar
%     set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
%     axis equal, axis xy
%     % abv = abs(max(units(idx(iunit)).smthRF(:)));
%     % colormap(redblue), clim([-abv abv])
%     pause
%     close
% end
% %

%% plot multiple at a time
%
doplot=true;
if doplot
    % idx=1:numel(units);
    idx = find([units.firing_rate]>2);% & [units.shank]==3);

    tic
    for iunit = 1:numel(idx)
        % For every new set of panels, open a new figure (except for the first set)
        if mod(iunit-1, 35) == 0
            % if iunit > 1
            % pause;      % Wait for user input after a set of 20
            % close(gcf); % Close the current figure
            % end
            figure;         % Open a new figure
        end

        % Create subplot in a 4x5 grid
        subplot(7, 5, mod(iunit-1, 35) + 1)

        % Plot the receptive field
        thisMap_s = units(idx(iunit)).smthRF;
        imagesc(flipud(thisMap_s))
        xlabel('Azimuth')
        ylabel('Elevation')
        xticks(linspace(1, size(thisMap_s, 2), 13))
        xticklabels(-120:20:120)
        yticks(linspace(1, size(thisMap_s, 1), 7))
        yticklabels(-60:20:60)
        title(sprintf('Unit %i', idx(iunit)))
        colorbar
        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 8)
        % axis equal,
        axis xy
        abv = abs(max(units(idx(iunit)).smthRF(:)));
        % colormap(redblue), clim([-abv abv])
    end

    toc
end

%% plot locations of good rfs

good_idx = [41 108, 145, 134, 156, 165, 175, 188, 283, 362, 367, 375, 498, 508, 600, 673, 752]
%
allocs= cat(1,units.location);

figure, hold on
plot(allocs(:,1), allocs(:,2),'k.');
plot(allocs([good_idx],1), allocs([good_idx],2),'r*')


pitch = 32;
1500/32



%% look at a subset of units


sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 100; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = 1:12;  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 1; % 1 if you want to see details of the receptive fields

bin_size = 1/100;
AnalysisTimeWindow = [0 1/100*12];


tic
for iunit = 1:length(units)
    spikes_this_cell = units(iunit).spike_times;

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, stim_times, AnalysisTimeWindow, bin_size);
    resps(iunit,:,:) = binnedArray';

end
toc

% for iunit = 1:numel(idx)
%     SUA_data = squeeze(resps(idx(iunit),:,:))'; % returns 1 x N time bins x nFrames;
% 
%     sn_options.figName =  ['sh: ', num2str(round(locs(idx(iunit),1)./250))];
%     tempRF = sparseNoiseAnalysis(stim_matrix,SUA_data,[],[],sn_options);
%     pause
% end
