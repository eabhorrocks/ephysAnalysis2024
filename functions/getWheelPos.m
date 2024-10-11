function wheel = getWheelPos(nidqDir,wheelChans)

nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

%% Read Wheel Data %%
wheel = NidqBin(wheelChans,:);

Wdata1=wheel(1,:);
Wdata2=wheel(2,:);

thresh1 = (max(Wdata1)-min(Wdata1))/2;
thresh2 = (max(Wdata2)-min(Wdata2))/2;
Wdata1_diff = diff(Wdata1); Wdata1_diff = [Wdata1_diff(1),Wdata1_diff];
Wdata2_diff = diff(Wdata2); Wdata2_diff = [Wdata2_diff(1),Wdata2_diff];

% Find rising edges
plocs_1 = find(Wdata1_diff > thresh1);
plocs_2 = find(Wdata2_diff > thresh2);
% Find falling edges
nlocs_1 = find(Wdata1_diff < -thresh1);
nlocs_2 = find(Wdata2_diff < -thresh2);

% this is based on arduino code -> not sure correct quadrature encoding?
temp_Position = zeros(size(Wdata1));
temp_Position(plocs_1(Wdata2(plocs_1) < thresh2)) = 1; %1
temp_Position(plocs_1(Wdata2(plocs_1) > thresh2)) = 0; %-1
temp_Position(nlocs_1(Wdata2(nlocs_1) < thresh2)) = 0; %-1
temp_Position(nlocs_1(Wdata2(nlocs_1) > thresh2)) = 0; %1

temp_Position(plocs_2(Wdata1(plocs_2) < thresh1)) = -1; %-1
temp_Position(plocs_2(Wdata1(plocs_2) > thresh1)) = 0; %1
temp_Position(nlocs_2(Wdata1(nlocs_2) < thresh1)) = 0; %1
temp_Position(nlocs_2(Wdata1(nlocs_2) > thresh1)) = 0; %-1

% Convert rotations to speed
wheel=struct;
wheel.pos = cumsum(temp_Position); %
wheel.time = (1:length(wheel.pos))./str2double(NidqMeta.niSampRate);


wheel.pos=wheel.pos(:);
wheel.time=wheel.time(:);

