function [t_ref_clean, t_sig_clean] = clean_pulse_vectors(t_pulse_ref, t_pulse_sig)
% CLEAN_PULSE_VECTORS Identifies and removes extra pulses from the longer vector.
%
% Inputs:
%   t_pulse_ref - Reference time vector
%   t_pulse_sig - Signal time vector
%
% Outputs:
%   t_ref_clean - The reference vector with extras removed
%   t_sig_clean - The signal vector (unchanged, unless it was the long one)

    % 1. Identify which is Long vs Short
    % We need to track which input is which so we can assign them back correctly at the end.
    if length(t_pulse_ref) > length(t_pulse_sig)
        t_long = t_pulse_ref;
        t_short = t_pulse_sig;
        long_is_ref = true; % Flag to remember Ref is the long one
    else
        t_long = t_pulse_sig;
        t_short = t_pulse_ref;
        long_is_ref = false; % Flag to remember Sig is the long one
    end

    % 2. Determine a safe threshold
    if length(t_long) > 1
        dt_avg = median(diff(t_long));
        threshold = dt_avg / 2;
    else
        threshold = 0.001; % Fallback
    end

    % 3. Find the extra pulse(s)
    extra_indices = [];
    
    for i = 1:length(t_long)
        % Find distance to the closest point in the short vector
        min_dist = min(abs(t_short - t_long(i)));
        
        % If the closest point is too far, it's an extra pulse
        if min_dist > threshold
            extra_indices = [extra_indices, i];
        end
    end
    
    % 4. Clean the vector (Remove the extras)
    if ~isempty(extra_indices)
        t_long(extra_indices) = []; % Deletes the elements at these indices
        fprintf('Removed %d pulse(s) from the longer vector.\n', length(extra_indices));
    else
        fprintf('No extra pulses found to remove.\n');
    end

    % 5. Assign back to the correct output variables
    if long_is_ref
        t_ref_clean = t_long;  % This is the cleaned version
        t_sig_clean = t_short; % This one was untouched
    else
        t_ref_clean = t_short; % This one was untouched
        t_sig_clean = t_long;  % This is the cleaned version
    end

end