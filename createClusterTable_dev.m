kilosort_directory = 'Z:\ibn-vision\DATA\SUBJECTS\M24019\ephys\20240719\spike_sorting\probe0\sorters\kilosort4\sorter_output';

% To do:
% check if any files are missing
% read sample_rate from meta file.

meta_file = dir(fullfile(kilosort_directory,'*tcat.imec0.ap.meta'));

st_file = dir(fullfile(kilosort_directory,['spike_times.npy']));
clusterID_file = dir(fullfile(kilosort_directory,['spike_clusters.npy']));
chanPos_file = dir(fullfile(kilosort_directory,['channel_positions.npy']));
% clusterTable_file = dir(fullfile(kilosort_directory,['clus_Table*.npy']));
% clusterTable_file = clusterTable_file(end); % get most recent

clusterNoise_file = dir(fullfile(kilosort_directory,['cluster_group.tsv'])); % noise not noise
clusterKSlabel_file = dir(fullfile(kilosort_directory,['cluster_KSLabel.tsv']));


%%

cluster_table = table;


sample_rate = 30000; % change to read from .meta file
% load spike times in samples and convert to seconds
spike_times_samples = double(readNPY(fullfile(kilosort_directory, st_file.name)));
spike_times = spike_times_samples/sample_rate;



% get cluster_id for each spike time
cluster_id = double(readNPY(fullfile(kilosort_directory,  clusterID_file.name)));
uniqueClusters = unique(cluster_id);

% add table var for spike times


cluster_table = addvars(cluster_table, uniqueClusters, 'NewVariableNames', 'cluster_id');
% unique cluster ids in the cluster_table

for icluster = 1:numel(uniqueClusters)
    table_row = find(cluster_table.cluster_id==uniqueClusters(icluster)); % what is the table row for this cluster_id?
    tempSpikeTimes = spike_times(cluster_id==uniqueClusters(icluster)); % get spike_times that correspond to this cluster_id
    cluster_table.spike_times(table_row) = {tempSpikeTimes}; % add these spike_times to the correct table_row
end


%%
% allenMetrics_file = dir(fullfile(kilosort_directory,['metrics*.csv']));
% allenMetrics_file = allenMetrics_file(end); % get most recent.
% 
% allenWaveforms_file = dir(fullfile(kilosort_directory,['waveform_metrics*.csv']));
% allenWaveforms_file = allenWaveforms_file(end);
% 
% meanWaveforms_file = dir(fullfile(kilosort_directory,['mean_waveforms*.npy']));
% meanWaveforms_file = meanWaveforms_file(end); % get most recent