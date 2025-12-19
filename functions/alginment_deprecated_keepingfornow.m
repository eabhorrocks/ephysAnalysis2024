%% 1. Preparation
pd_onsets = nidq.imecTime(final_on_indices);
arduino_onsets = trialParams.ArduinoTime / 1000; % Convert to seconds

% Calculate Inter-Onset Intervals (IOI) - The "Fingerprint" of your session
pd_IOI = diff(pd_onsets); 
pd_IOI = pd_IOI(:); % Ensure column vector

arduino_IOI = diff(arduino_onsets);
arduino_IOI = arduino_IOI(:); % Ensure column vector

%% 2. Diagnose the Mismatch
n_pd = length(pd_IOI);
n_arduino = length(arduino_IOI);

fprintf('Photodiode Intervals: %d\n', n_pd);
fprintf('Arduino Intervals:    %d\n', n_arduino);

if n_pd > n_arduino
    fprintf('WARNING: Found %d EXTRA pulses in the photodiode signal.\n', n_pd - n_arduino);
    
    %% 3. Find the Best Alignment (Cross-Correlation)
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
    
  %% Apply the Alignment (Offset 1)
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

%% Assign to Table
% We use (:) to force them into Column vectors, preventing dimension errors
trialParams.PD_OnsetIndex = final_pd_indices(:);
trialParams.PD_OffsetIndex = final_pd_offsets(:);
trialParams.PD_OnsetVal = final_pd_times(:);
trialParams.PD_ActualDuration = final_pd_durations(:);

fprintf('Success: Assigned %d aligned photodiode trials to the table.\n', height(trialParams));

%% Optional: Final Sanity Check Plot
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