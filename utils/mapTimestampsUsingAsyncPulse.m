function [timestamps_mapped] = mapTimestampsUsingAsyncPulse(syncTimes_input,syncTimes_mapto,timestamps_input)

% times must be in the same unit

%% use xcorr to find matching sync pulses from the two data streams
[r, lags] = xcorr(diff(syncTimes_input), diff(syncTimes_mapto));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);
% ...then remove the syncpulses from bonsai trace that occur before the
% optimum lag
if best_lag > 0
    syncTimes_input = syncTimes_input(syncTimes_input>syncTimes_input(best_lag));
end

% ...and rerun the xcorr
[r, lags] = xcorr(diff(syncTimes_input), diff(syncTimes_mapto));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);

% check which async pulse sequence is available first (ie. which was turned
% on first)
if best_lag < 0
    nSyncOffset = -best_lag+1;
    syncTimes_mapto_processed = syncTimes_mapto(nSyncOffset:end); % sync sglx times
    if length(syncTimes_mapto_processed) > length(syncTimes_input)
        syncTimes_mapto_processed = syncTimes_mapto_processed(1:length(syncTimes_input));
    end
    syncTimes_input_processed = syncTimes_input(1:numel(syncTimes_mapto_processed)); %sync bonsai times
else
    nSyncOffset = best_lag+1;
    syncTimes_mapto_processed = syncTimes_mapto(1:end);
    if length(syncTimes_mapto_processed) > length(syncTimes_input)
        syncTimes_mapto_processed = syncTimes_mapto_processed(1:length(syncTimes_input));
    end
    syncTimes_input_processed = syncTimes_input(nSyncOffset:numel(syncTimes_mapto_processed)+nSyncOffset-1);
end

%% use aligned sync pulses to map timestamps

timestamps_mapped = interp1(syncTimes_input_processed, syncTimes_mapto_processed, timestamps_input,'linear','extrap');

end
