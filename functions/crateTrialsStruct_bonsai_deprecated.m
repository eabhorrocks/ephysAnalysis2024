function trials = crateTrialsStruct_bonsai_deprecated(nidq, csvFilePath)
%crateTrialsStruct_flicker Aligns trial data and calculates flicker phase vectors at 1 kHz.
%
%   This function intelligently aligns trials based on stimulus frequency, making it robust
%   to missing trials. It calculates the continuous phase vector in radians
%   for two color channels at a resampled rate of 1 kHz.
%
%   SYNTAX:
%   trials = crateTrialsStruct_flicker(nidq, csvFilePath)
%
%   INPUTS:
%   nidq         - A struct from extractNidqData, containing 'stimON', 'stimCHANGE', etc.
%   csvFilePath  - Path to the CSV file with trial parameters. May have a header.
%
%   OUTPUT:
%   trials       - A struct array where each element represents an aligned trial,
%                  containing parameters, measured data, and phase vectors.

%% 1. Load Data Sources
% Use imecTime if available, otherwise default to NIDQ time
if isfield(nidq, 'imecTime') && ~isempty(nidq.imecTime)
    timeVector = nidq.imecTime;
    timeBaseSource = 'imecTime';
else
    timeVector = nidq.time;
    timeBaseSource = 'NIDQ time';
end

% Load trial parameters from CSV, automatically handling a potential header
if ~isfile(csvFilePath)
    error('CSV file not found at: %s', csvFilePath);
end
opts = detectImportOptions(csvFilePath);
% Standardize variable names to use 'green' and 'uv'
% opts.VariableNames = {'stimType', 'duration_ms', 'frequency_param', 'phase_green', 'phase_uv', ...
%     'contrast_green', 'contrast_uv', 'bonsaiTime', 'arduinoTime'};
trialParams = readtable(csvFilePath, opts);
num_csv_trials = height(trialParams);
fprintf('Loaded %d trials from parameter file.\n', num_csv_trials);

%% 2. Pre-process NIDQ Data to Find All Potential Trials

% Photodiode signal can occasionally glitch creating false LOW in the
% middle of a stimulus, so we want to detect and fill these in
% appropriately.
d_sig = diff([0, nidq.photodiode, 0]); % Pad to handle edges safely
on_indices  = find(d_sig == 1);
off_indices = find(d_sig == -1) - 1; % -1 to correct index shift from padding

% Clean up incomplete trials at the edges
if ~isempty(off_indices) && ~isempty(on_indices)
    if off_indices(1) < on_indices(1), off_indices(1) = []; end
    if ~isempty(on_indices) && on_indices(end) > off_indices(end), on_indices(end) = []; end
end

% Safety check
n_pulses = min(length(on_indices), length(off_indices));
on_indices = on_indices(1:n_pulses);
off_indices = off_indices(1:n_pulses);

% Identify and Fill Gaps (transient LOWs)
min_gap_sec = 0.01; % 10ms threshold

% Get duration of gaps in time using your imecTime vector
t_off = nidq.imecTime(off_indices(1:end-1));
t_on_next = nidq.imecTime(on_indices(2:end));
gaps_sec = t_on_next - t_off;

% Find which gaps are "glitches" (shorter than min_gap)
glitch_indices = find(gaps_sec < min_gap_sec);

fprintf('Found %d transient drops (< %.3fs) to bridge.\n', length(glitch_indices), min_gap_sec);

% Loop through and fill them with 1s
% We iterate backwards to avoid any indexing confusion, though not strictly necessary here
for i = length(glitch_indices):-1:1
    idx = glitch_indices(i);
    
    % The gap is between off_indices(idx) and on_indices(idx+1)
    % We fill from the offset to the next onset
    fill_start = off_indices(idx); 
    fill_end = on_indices(idx+1);
    
    nidq.photodiode(fill_start:fill_end) = 1;
end

% Identify transient highs and fill with 0s
% We need to compare the current ON times to the current OFF times.
t_on_current = nidq.imecTime(final_on_indices);
t_off_current = nidq.imecTime(final_off_indices);

pulse_durations = t_off_current - t_on_current;

% Find which pulses are "spikes" (shorter than min_gap_sec)
spike_indices = find(pulse_durations < min_gap_sec);
fprintf('Found %d transient spikes (< %.3fs) to remove.\n', length(spike_indices), min_gap_sec);

% Loop through and set them to 0 (remove the pulse)
for i = 1:length(spike_indices)
    idx = spike_indices(i);
    
    % Get the start/end samples of this specific short pulse
    clear_start = final_on_indices(idx);
    clear_end = final_off_indices(idx);
    
    nidq.photodiode(clear_start:clear_end) = 0;
end




% 3. Final Edge Detection (on CLEAN signal)
% Now that the gaps are filled, the split trials have merged into single long trials.
% We re-run the detection to get the final, clean timestamps.

d_sig_clean = diff([0, nidq.photodiode, 0]);
final_on_indices = find(d_sig_clean == 1);
final_off_indices = find(d_sig_clean == -1) - 1;

% Clean up edges one last time (standard safety)
if ~isempty(final_off_indices) && ~isempty(final_on_indices)
    if final_off_indices(1) < final_on_indices(1)
        final_off_indices(1) = []; 
    end
    if ~isempty(final_on_indices) && final_on_indices(end) > final_off_indices(end)
        final_on_indices(end) = []; 
    end
end

% now we want to check whether some pause (or other) trials were presented
% in bonsai after the ephys stopped. To do so we check the trials
% beforehand are aligned in their timing, and then remove the excess trials
% from the trialParams table.
pd_onsets = nidq.imecTime(final_on_indices);
num_detected = length(pd_onsets);
arduino_onsets = trialParams.ArduinoTime/1000; % convert ms to s

pd_IOI = diff(pd_onsets); pd_IOI=pd_IOI(:);
arduino_IOI = diff(arduino_onsets);

n_pd = length(pd_IOI);
n_arduino = length(arduino_IOI);

fprintf('Photodiode Intervals: %d\n', n_pd);
fprintf('Arduino Intervals:    %d\n', n_arduino);

arduino_IOI_truncated = arduino_IOI(1:numel(pd_IOI));
time_diff = abs(pd_IOI - arduino_IOI_truncated);

tolerance_vec = repmat(0.1, size(time_diff)); 
tolerance_vec(1) = 0.5; % Set the first trial's tolerance to the lenient value (500ms)

is_aligned_per_trial = time_diff < tolerance_vec;
is_aligned_overall = all(is_aligned_per_trial);
mean_error = mean(time_diff);


% cross-corr alignment

if n_pd > n_arduino
    fprintf('WARNING: Found %d EXTRA pulses in the photodiode signal.\n', n_pd - n_arduino);
    
    % Find the Best Alignment (Cross-Correlation)
    % We assume the Arduino sequence is the "template" that exists somewhere 
    % inside the noisy Photodiode sequence.
    
    % We sweep the Arduino sequence across the PD sequence to find the match.
    % To do this efficiently with xcorr, we center the data.
    
    % Note: xcorr can be messy with varying magnitudes. 
    % A direct "Sum of Squared Differences" (SSD) search is often cleaner for timing.
    
    min_error = inf;
    best_offset = 0;
    
    % We only need to check offsets that allow the Arduino seq to fit inside PD seq
    max_offset = n_pd - n_arduino;
    
    for offset = 0:max_offset
        % Extract the snippet from PD that matches the length of Arduino
        pd_snippet = pd_IOI( (1:n_arduino) + offset );
        
        % Calculate error (L1 norm: sum of absolute differences)
        current_error = sum(abs(pd_snippet - arduino_IOI));
        
        if current_error < min_error
            min_error = current_error;
            best_offset = offset;
        end
    end
    
    fprintf('Best fit found at index offset: %d\n', best_offset);
    fprintf('Total timing error at this alignment: %.4f seconds\n', min_error);
    
% Apply the Alignment (Offset 1)
% Recalculate the indices using the found offset
n_trials_table = height(trialParams); % e.g., 1618
valid_start_idx = best_offset + 1;    % Start at pulse 2
valid_end_idx = best_offset + n_trials_table; % Grab next 1618 pulses

% Extract the correct range of indices
final_pd_indices = final_on_indices(valid_start_idx : valid_end_idx);
final_pd_offsets = final_off_indices(valid_start_idx : valid_end_idx);

% Get the actual timestamps
final_pd_times = nidq.imecTime(final_pd_indices);
final_pd_durations = nidq.imecTime(final_pd_offsets) - final_pd_times;

% Assign to Table
% We use (:) to force them into Column vectors, preventing dimension errors
trialParams.PD_OnsetIndex = final_pd_indices(:);
trialParams.PD_OffsetIndex = final_pd_offsets(:);
trialParams.PD_OnsetVal = final_pd_times(:);
trialParams.PD_ActualDuration = final_pd_durations(:);

fprintf('Success: Assigned %d aligned photodiode trials to the table.\n', height(trialParams));

% Optional: Final Sanity Check Plot
% Plot the duration mismatch. If alignment is wrong, this will look like garbage.
% If alignment is right, this will show your true display variability.
figure; 
plot(trialParams.Duration - trialParams.PD_ActualDuration, '.');
yline(0, 'r');
ylabel('Error (s)'); xlabel('Trial #');
title('Duration Error (Expected - Actual)');
% If you see outliers here, those are real specific trials where the frame dropped.
    
else
    disp('Problem is not "more PD than Arduino". Check previous steps.');
end



% optional visual check for debugging
% figure('Name', 'Trial Alignment Check');
% subplot(2,1,1);
% plot(normalize(pd_IOI), 'b', 'LineWidth', 1.5, 'DisplayName', 'NIDAQ');
% hold on;
% plot(normalize(arduino_IOI_truncated), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Bonsai');
% title('Rhythm Comparison');
% legend; grid on;
% 
% subplot(2,1,2);
% plot(time_diff, 'k', 'DisplayName', 'Actual Error');
% hold on;
% % Plot the tolerance threshold in red to visualize pass/fail
% plot(tolerance_vec, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Allowed Tolerance');
% title('Timing Discrepancy vs Tolerance');
% ylabel('Error (s)'); xlabel('Trial Interval Index');
% legend; grid on;

if is_aligned_overall
    fprintf('\nSUCCESS: Pulses match (accounting for startup jitter).\n');
    fprintf('Truncating last %d rows from trialParams.\n', height(trialParams) - num_detected);
    
    % Truncate table
    trialParams = trialParams(1:num_detected, :);
    
    % Assign precise timings
    start_times = nidq.imecTime(final_on_indices(:));
    end_times = nidq.imecTime(final_off_indices(:));
    trialParams.startTime = start_times(:);
    trialParams.endTime = end_times(:);
else
    % If it fails, identify WHERE it failed
    first_fail = find(~is_aligned_per_trial, 1);
    warning('FAILURE: Rhythm mismatch at Interval #%d (Error: %.4fs)', first_fail, time_diff(first_fail));
end

trials = table2struct(trialParams);




%% helper function


