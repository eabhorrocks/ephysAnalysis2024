function trials = alignDataToTrials(trials, dataTimeVector, preStim_s, postStim_s, varargin)
%alignDataToTrials Aligns continuous data streams to trial event times.
%
%   This function takes a trial structure and continuous data (e.g., wheel
%   speed, eye position) and extracts snippets for each trial corresponding
%   to a window around the stimulus. It adds these snippets and their mean
%   values as new fields to the input trial structure.
%
%   SYNTAX:
%   trials = alignDataToTrials(trials, dataTimeVector, preStim_s, postStim_s, ...
%                              'dataName1', dataVector1, 'dataName2', dataVector2, ...)
%
%   INPUTS:
%   trials         - A struct array where each element must contain at least
%                    a 'startTime' and 'endTime' field in seconds.
%   dataTimeVector - The master time vector (in seconds) for the continuous data.
%   preStim_s      - The duration (in seconds) to include before stimulus onset.
%   postStim_s     - The duration (in seconds) to include after stimulus offset.
%   varargin       - Name-value pairs of the data to align. The 'name'
%                    (a string) will be used to name the new fields in the
%                    output struct. The 'value' is the corresponding data vector,
%                    which must be the same length as dataTimeVector.
%
%   OUTPUT:
%   trials         - The input struct array with new fields added for each
%                    data stream. For each input 'dataName', the following
%                    fields are added to each trial:
%                       .time_aligned: Time vector for the snippet, zeroed to stim onset.
%                       .dataName: The data snippet from pre- to post-stim.
%                       .dataName_mean: The mean of the full data snippet (NaNs ignored).
%                       .dataName_stim_period: The data snippet during the stim period only.
%                       .dataName_stim_period_mean: The mean of the stim-period data (NaNs ignored).

%% 1. Input Parsing and Validation
if mod(length(varargin), 2) ~= 0
    error('Data inputs must be provided as name-value pairs.');
end

numDataStreams = length(varargin) / 2;
dataStreams = struct();
for i = 1:numDataStreams
    dataName = varargin{2*i - 1};
    dataVector = varargin{2*i};

    if ~ischar(dataName)
        error('Input data names must be character strings.');
    end
    if length(dataVector) ~= length(dataTimeVector)
        error('Data vector ''%s'' must be the same length as the dataTimeVector.', dataName);
    end
    
    % Store the validated data and name
    dataStreams(i).Name = dataName;
    dataStreams(i).Vector = dataVector;
end

fprintf('Aligning %d data stream(s) to %d trials...\n', numDataStreams, length(trials));

%% 2. Process Each Trial
for i = 1:length(trials)
    % Skip trials that were not successfully matched and lack timing info
    if isnan(trials(i).startTime)
        % Add empty fields for consistency
        trials(i).time_aligned = [];
        for j = 1:numDataStreams
            dataName = dataStreams(j).Name;
            trials(i).(dataName) = [];
            trials(i).([dataName '_mean']) = NaN; % Add mean field as NaN
            trials(i).([dataName '_stim_period']) = [];
            trials(i).([dataName '_stim_period_mean']) = NaN; % Add mean field as NaN
        end
        continue;
    end

    % --- Define Time Windows for this Trial ---
    stim_onset = trials(i).startTime;
    stim_offset = trials(i).endTime;
    
    % Window for the full snippet (pre-stimulus to post-stimulus)
    snippet_start_time = stim_onset - preStim_s;
    snippet_end_time = stim_offset + postStim_s;

    % --- Extract Full Snippet ---
    snippet_indices = dataTimeVector >= snippet_start_time & dataTimeVector <= snippet_end_time;
    
    % Create and add the zero-aligned time vector for this snippet
    trials(i).time_aligned = dataTimeVector(snippet_indices) - stim_onset;

    % --- Extract Stimulus-Period-Only Snippet ---
    stim_period_indices = dataTimeVector >= stim_onset & dataTimeVector <= stim_offset;

    % --- Loop Through Each Data Stream and Add to Struct ---
    for j = 1:numDataStreams
        dataName = dataStreams(j).Name;
        dataVector = dataStreams(j).Vector;
        
        % --- Full Snippet ---
        % Extract data
        snippet_data = dataVector(snippet_indices);
        % Store data snippet
        trials(i).(dataName) = snippet_data;
        % Calculate and store mean (ignoring NaNs)
        trials(i).([dataName '_mean']) = mean(snippet_data, 'omitnan');

        % --- Stimulus Period Only ---
        % Extract data
        stim_period_data = dataVector(stim_period_indices);
        % Store data snippet
        trials(i).([dataName '_stim_period']) = stim_period_data;
        % Calculate and store mean (ignoring NaNs)
        trials(i).([dataName '_stim_period_mean']) = mean(stim_period_data, 'omitnan');
    end
end

fprintf('Alignment complete.\n');

end