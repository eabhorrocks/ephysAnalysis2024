%% create cluster table (DEV)

function units = createClusterTable(ephysDir, SorterName)

% Subject = 'M24019';
% Session = '20240718';
% BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% SorterName = 'kilosort4_merged';
kilosortDir = fullfile(ephysDir,'spike_sorting','probe0','sorters',SorterName);
extensionsDir = fullfile(ephysDir,'spike_sorting','probe0','sorters',SorterName,'waveform','extensions');


imSampRate=30000; % load from .meta file


%% get pre-computed quality metrics
quality_metrics = readtable(fullfile(extensionsDir,'quality_metrics','metrics.csv'));
template_metrics = readtable(fullfile(extensionsDir,'template_metrics','metrics.csv'));
cluster_metrics = join(quality_metrics,template_metrics);
cluster_metrics = renamevars(cluster_metrics,"Var1","cluster_id");


%% get spike times

spikes = readmatrix(fullfile(kilosortDir,'spikes.csv'));
spike_samples = spikes(:,3);
spike_times = spike_samples./imSampRate;
spike_clusters = spikes(:,1);
original_cluster_id = unique(spike_clusters);
spike_clusters_updated = cluster_metrics.cluster_id(spike_clusters+1);

for icluster = 1:height(cluster_metrics)
    cluster_metrics.spike_times(icluster) = {spike_times(spike_clusters_updated==cluster_metrics.cluster_id(icluster))};
end


%% get waveforms and peak channel

templates = readNPY(fullfile(extensionsDir,'templates','average.npy'));
[~,peak_channel] = max(squeeze(mean(abs(templates),2)),[],2); % channel index providing the maximum of the mean of the absolute
cluster_metrics.peak_channel = peak_channel;

for icluster = 1:height(cluster_metrics)
    peak_channel_waveforms(icluster,:) = templates(icluster,:,peak_channel(icluster));
    cluster_metrics.peak_channel_waveform(icluster,:) = peak_channel_waveforms(icluster,:);
end

%% get unit locations in co-ordinates

unitLocs = readNPY(fullfile(extensionsDir,'unit_locations','unit_locations.npy'));
for icluster = 1:height(cluster_metrics)
    cluster_metrics.location(icluster,:) = unitLocs(icluster,:);
end

%%

units = cluster_metrics;
% units = table2struct(cluster_metrics);

end
%% get physiologically-defined cell type
% 
% tic
% 
% for icluster = 1:height(cluster_metrics)
%     [ccg, t] = CCG(cluster_metrics.spike_times{icluster}, ones(size(cluster_metrics.spike_times{icluster})),...
%     'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
%     cluster_metrics.acg(icluster,:) = ccg';
%     fit_params_out = fit_ACG(ccg,false);
%     cluster_metrics.tau_rise(icluster) = fit_params_out.acg_tau_rise;
% end
% toc
