function sta = compute_sta_multi_frame_avg_spikes(stimulus, frame_times, spike_times, nLags)
% compute_sta_multi_frame_avg_spikes computes the multi-frame spike-triggered average (STA)
% by first extracting, for each spike, the preceding nLags frames and then averaging these 
% individual STAs over spikes.
%
% Inputs:
%   stimulus    - 3D matrix of size [nAzi, nEle, nFrames]. The first two dimensions 
%                 are spatial (e.g., x and y pixels) and the third indexes frames.
%   frame_times - Vector of length nFrames containing the time at which each frame was shown.
%   spike_times - Vector containing the spike times.
%   nLags       - Number of frames to include before each spike (e.g., 10).
%
% Output:
%   sta         - A 3D matrix of size [nAzi, nEle, nLags] representing the averaged 
%                 STA. For example, sta(:,:,1) is the average stimulus frame immediately 
%                 preceding the spike, sta(:,:,2) is the average for the frame one step earlier, etc.

[nAzi, nEle, nFrames] = size(stimulus);
sta_sum = zeros(nAzi, nEle, nLags); % Accumulator for STA from each spike
valid_spike_count = 0;              % Counter for spikes with enough preceding frames

for i = 1:length(spike_times)
    spike = spike_times(i);
    % Find the index of the last frame presented before the spike.
    idx = find(frame_times <= spike, 1, 'last');
    
    % Only consider spikes that have at least nLags preceding frames.
    if ~isempty(idx) && (idx >= nLags)
        % Extract the nLags frames preceding the spike.
        spike_sta = stimulus(:, :, idx - nLags + 1 : idx);
        % Add the spike's STA to the accumulator.
        sta_sum = sta_sum + spike_sta;
        valid_spike_count = valid_spike_count + 1;
    end
end

if valid_spike_count > 0
    % Average the accumulated STAs over the valid spikes.
    sta = sta_sum / valid_spike_count;
else
    warning('No spikes had enough preceding frames to compute the STA.');
    sta = [];
end
end
