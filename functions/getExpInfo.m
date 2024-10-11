function [ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session)

% function uses a gkey file (user-made) to construct an Experiment Info
% table with corresponding sample/timestamps and bonsai csv files for each
% stimulus

%% paths to required directories
ephysDir = fullfile(BaseDir,Subject,'ephys',Session);
nidqDir = fullfile(BaseDir,Subject,'ephys',Session,'nidq_processed');
bonsaiDir = fullfile(BaseDir,Subject,'Bonsai',Session);

%% load gkey + catGT files and generate exp info table

gkeyFile = dir(fullfile(ephysDir,'*gkey.csv'));
gkeyTable = readtable(fullfile(gkeyFile.folder, gkeyFile.name));

% sample/time offsets for catGT concatenation
ctOffsetFiles = dir(fullfile(nidqDir,'*ct_offsets.txt'));
scOffsetFiles = dir(fullfile(nidqDir,'*sc_offsets.txt'));

% check if sueprcat done, otherwise create dummy table
if numel(scOffsetFiles)>=1
    s = rows2vars(readtable(fullfile(scOffsetFiles(1).folder,scOffsetFiles(1).name)));
    s(:,1) = [];
    for ivar = 1:width(s)
        varNames(ivar) = s{1,ivar};
    end
    s.Properties.VariableNames = varNames;
    s(1,:)=[];
    s = convertvars(s,s.Properties.VariableNames, 'cell2mat');
else
    s= table(0,0,0,0,0,0,'VariableNames',{'smp_imap0:', 'smp_nidq:', 'sec_imap0:', 'sec_nidq:', 'srate_imap0:', 'srate_nidq:'});
end
supercat_tbl = s; clear s

% loop through cat GT offset files
for ifile = 1:numel(ctOffsetFiles)
    s = rows2vars(readtable(fullfile(ctOffsetFiles(ifile).folder,ctOffsetFiles(ifile).name)));
    s(:,1) = [];

    varNames={};
    for ivar = 1:width(s)
        varNames(ivar) = s{1,ivar};
    end
    s.Properties.VariableNames = varNames;
    s(1,:)=[];
    s = convertvars(s,s.Properties.VariableNames, 'cell2mat');
    for ivar = 1:width(s)
        s.(varNames{ivar}) = s.(varNames{ivar})+supercat_tbl.(varNames{ivar})(ifile);
    end

    % get acq name
    temp = strsplit(ctOffsetFiles(ifile).name,'_');
    acqName = [temp{1},'_',temp{2},'_',temp{3}];
    acqCol = repmat(acqName,height(s),1);
    s.Acquisition = acqCol;

    tempTbl{ifile} = s; clear s
end

tempStimInfo = cat(1,tempTbl{:});

% create stim info table, combining above, with nice ordering
ExpInfo = table;
ExpInfo.AcquisitionName = cellstr(tempStimInfo.Acquisition);
ExpInfo.AcqusitionIndex = gkeyTable.Acquisiton_index;
ExpInfo.gKey = gkeyTable.gNum;
ExpInfo.smp_imap0 = tempStimInfo.('smp_imap0:');
ExpInfo.sec_imap0 = tempStimInfo.('sec_imap0:');
ExpInfo.smp_nidq = tempStimInfo.('smp_nidq:');
ExpInfo.sec_nidq = tempStimInfo.('sec_nidq:');
ExpInfo.Stimulus = gkeyTable.Stimulus;
ExpInfo.notes = gkeyTable.notes;

emptyCol = cell(height(ExpInfo),1);
ExpInfo.TrialParamsFiles = emptyCol;
ExpInfo.WheelFiles = emptyCol;
ExpInfo.PhotodiodeFiles = emptyCol;
ExpInfo.VideoCSVFiles = emptyCol;
ExpInfo.VideoFiles = emptyCol;


%% find unique stimuli shown and get corresponding bonsai files

uniqueStimuli = unique(ExpInfo.Stimulus);

for istim = 1:numel(uniqueStimuli)

    stimulus = uniqueStimuli{istim};

    trialParamsFiles = dir(fullfile(bonsaiDir,['*',stimulus,'_TrialParams*']));
    wheelFiles = dir(fullfile(bonsaiDir,['*',stimulus,'_Wheel*']));
    photodiodeFiles = dir(fullfile(bonsaiDir,['*',stimulus,'_Photodiode*']));
    videoCSVFiles = dir(fullfile(bonsaiDir,['*',stimulus,'_video*.csv']));
    videoFiles = dir(fullfile(bonsaiDir,['*',stimulus,'_video*.avi']));

    nParamsFiles = numel(trialParamsFiles);
    nWheelFiles = numel(trialParamsFiles);
    nPhotodiodeFiles = numel(trialParamsFiles);
    nVideoCSVFiles = numel(trialParamsFiles);
    nVideoFiles = numel(trialParamsFiles);

    nGkeyIdx = sum(strcmp(ExpInfo.Stimulus, stimulus));

    if ~strcmp(stimulus,'Spontaneous')

    if nParamsFiles~=nGkeyIdx
        error(['Number of params files found doesnt match number of gkey entries for stimulus: ', stimulus])
    end
    if nWheelFiles~=nGkeyIdx
        error(['Number of wheel files found doesnt match number of gkey entries for stimulus: ', stimulus])
    end
    if nPhotodiodeFiles~=nGkeyIdx
        error(['Number of photodiode files found doesnt match number of gkey entries for stimulus: ', stimulus])
    end
    if nVideoCSVFiles~=nGkeyIdx
        error(['Number of video csv files found doesnt match number of gkey entries for stimulus: ', stimulus])
    end
    if nVideoFiles~=nGkeyIdx
        error(['Number of video files found doesnt match number of gkey entries for stimulus: ', stimulus])
    end

    allGkeyIdx = find((strcmp(ExpInfo.Stimulus, stimulus)));

    for ifile = 1:nGkeyIdx
        thisTableIndex = allGkeyIdx(ifile);

        ExpInfo.TrialParamsFiles(thisTableIndex) = {fullfile(trialParamsFiles(ifile).folder,trialParamsFiles(ifile).name)};
        ExpInfo.WheelFiles(thisTableIndex) = {fullfile(wheelFiles(ifile).folder,wheelFiles(ifile).name)};
        ExpInfo.PhotodiodeFiles(thisTableIndex) = {fullfile(photodiodeFiles(ifile).folder,photodiodeFiles(ifile).name)};
        ExpInfo.VideoCSVFiles(thisTableIndex) = {fullfile(videoCSVFiles(ifile).folder,videoCSVFiles(ifile).name)};
        ExpInfo.VideoFiles(thisTableIndex) = {fullfile(videoFiles(ifile).folder,videoFiles(ifile).name)};
    end

    end

end