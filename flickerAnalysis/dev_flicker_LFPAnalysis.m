%%
Subject = 'M24019';
Session = '20240718';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort4','sorter_output');
lfpDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0');
[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session);


%%
lfpDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0');
lfpMetaFile =  dir(fullfile(lfpDir,'*.meta'));
lfpBinFile =  dir(fullfile(lfpDir,'*.bin'));

lfpMeta = ReadMeta(fullfile(lfpFile.folder, lfpMetaFile.name));

dataArrayCorrect = ReadSingleChannelBin_edd(1, lfpMeta.imSampRate*10, lfpMeta, lfpBinFile.name, lfpBinFile.folder);

%%
channel = 377,327,189, 214;

dataArraySS = ReadBinEddit(channel, lfpMeta, lfpBinFile.name, lfpBinFile.folder)

%dataArray = ReadBinEddit(1, inf, channel, lfpMeta, lfpBinFile.name, lfpBinFile.folder)