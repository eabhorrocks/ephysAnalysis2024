function clusterTable = createClusterTable(kilosortDir)

% very dev until SI output sorted

%% get file paths
meta_file = dir(fullfile(kilosortDir,'*tcat.imec0.ap.meta'));

st_file = dir(fullfile(kilosortDir,['spike_times.npy']));
clusterID_file = dir(fullfile(kilosortDir,['spike_clusters.npy']));
chanPos_file = dir(fullfile(kilosortDir,['channel_positions.npy']));
% clusterTable_file = dir(fullfile(kilosortDir,['clus_Table*.npy']));
% clusterTable_file = clusterTable_file(end); % get most recent

clusterNoise_file = dir(fullfile(kilosortDir,['cluster_group.tsv'])); % noise not noise
clusterKSlabel_file = dir(fullfile(kilosortDir,['cluster_KSLabel.tsv']));

%% create cluster table

clusterTable = table;
sample_rate = 30000; % change to read from .meta file
% load spike times in samples and convert to seconds
spike_times_samples = double(readNPY(fullfile(kilosortDir, st_file.name)));
spike_times = spike_times_samples/sample_rate;

% get cluster_id for each spike time
cluster_id = double(readNPY(fullfile(kilosortDir,  clusterID_file.name)));
uniqueClusters = unique(cluster_id);

% add table var for spike times
clusterTable = addvars(clusterTable, uniqueClusters, 'NewVariableNames', 'cluster_id');
% unique cluster ids in the clusterTable

for icluster = 1:numel(uniqueClusters)
    table_row = find(clusterTable.cluster_id==uniqueClusters(icluster)); % what is the table row for this cluster_id?
    tempSpikeTimes = spike_times(cluster_id==uniqueClusters(icluster)); % get spike_times that correspond to this cluster_id
    clusterTable.spike_times(table_row) = {tempSpikeTimes}; % add these spike_times to the correct table_row
end