function [nidq_1hz_path, imec_1hz_path] = get_1hz_sync_files(nidqDir)
% GET_1HZ_SYNC_FILES  Locates or generates stitched 1Hz sync pulse files.
%
%   [nidqPath, imecPath] = GET_1HZ_SYNC_FILES(nidqDir)
%
%   INPUTS:
%       nidqDir - (String/Char) Full path to the directory containing 
%                 the .txt sync files.
%
%   OUTPUTS:
%       nidq_1hz_path - Full path to the valid NIDQ 1Hz text file.
%       imec_1hz_path - Full path to the valid IMEC 1Hz text file.
%
%   BEHAVIOR:
%       1. Scans nidqDir for standard sync files (*nidq.xd_8... and *imec0.ap...).
%       2. If 1 file pair found: Returns those paths directly.
%       3. If multiple file pairs found (SuperCat):
%          a. Checks if stitched files (*_sc_1hz.txt) already exist.
%          b. If not, loads *sc_offsets.txt, applies offsets, stitches, 
%             and saves the new files.
%          c. Returns paths to the stitched files.

    %% 1. Find Raw Files
    % Define patterns
    nidq_pattern = '*nidq.xd_8_3_500.txt';
    imec_pattern = '*imec0.ap.xd_384_6_500.txt';

    raw_nidq = dir(fullfile(nidqDir, nidq_pattern));
    raw_imec = dir(fullfile(nidqDir, imec_pattern));
    
    nFiles = numel(raw_nidq);

    % Sanity Checks
    if nFiles == 0
        error('get_1hz_sync_files:NoFiles', ...
            'No NIDQ files matching pattern "%s" found in %s', nidq_pattern, nidqDir);
    end
    
    if nFiles ~= numel(raw_imec)
        error('get_1hz_sync_files:Mismatch', ...
            'Mismatch: Found %d NIDQ files but %d IMEC files.', nFiles, numel(raw_imec));
    end

    %% 2. Process Files (Standard vs SuperCat)
    if nFiles == 1
        % --- Case A: Standard (Single File) ---
        fprintf('Standard Mode: Single file detected.\n');
        nidq_1hz_path = fullfile(raw_nidq.folder, raw_nidq.name);
        imec_1hz_path = fullfile(raw_imec.folder, raw_imec.name);
        
    else
        % --- Case B: SuperCat (Multiple Files) ---
        fprintf('SuperCat Mode: %d files detected.\n', nFiles);
        
        stitched_nidq_name = 'nidq_sc_1hz.txt';
        stitched_imec_name = 'imec_sc_1hz.txt';
        
        stitched_nidq_path = fullfile(nidqDir, stitched_nidq_name);
        stitched_imec_path = fullfile(nidqDir, stitched_imec_name);
        
        % Check if output already exists
        if exist(stitched_nidq_path, 'file') && exist(stitched_imec_path, 'file')
            fprintf('  > Stitched files found. Returning existing paths.\n');
            nidq_1hz_path = stitched_nidq_path;
            imec_1hz_path = stitched_imec_path;
        else
            % --- Stitching Logic ---
            fprintf('  > Stitched files missing. Generating...\n');
            
            % 1. Load Offsets
            offset_file = dir(fullfile(nidqDir, '*sc_offsets.txt'));
            if isempty(offset_file)
                error('get_1hz_sync_files:NoOffsets', 'SuperCat offsets file (*sc_offsets.txt) not found.');
            end
            
            sc_offsets = readtable(fullfile(offset_file.folder, offset_file.name));
            
            % 2. Extract Offsets Safely (Explicit Name Lookup)
            % Do not rely on row order. Explicitly find the row.
            nidq_row = strcmp(sc_offsets.Var1, 'sec_nidq:');
            imec_row = strcmp(sc_offsets.Var1, 'sec_imap0:');
            
            if ~any(nidq_row) || ~any(imec_row)
                 error('get_1hz_sync_files:BadOffsets', 'Could not find "sec_nidq:" or "sec_imap0:" in offsets table.');
            end
            
            % Extract columns 2:end
            nidq_offsets = table2array(sc_offsets(nidq_row, 2:end));
            imec_offsets = table2array(sc_offsets(imec_row, 2:end));
            
            % Check dimensions
            if numel(nidq_offsets) < nFiles || numel(imec_offsets) < nFiles
                 error('get_1hz_sync_files:OffsetMismatch', 'Offsets table has fewer columns than number of files.');
            end

            % 3. Stitch NIDQ
            fprintf('  > Stitching NIDQ... ');
            nidq_concat = [];
            for i = 1:nFiles
                fpath = fullfile(raw_nidq(i).folder, raw_nidq(i).name);
                raw_data = load(fpath);
                % Apply offset
                nidq_concat = [nidq_concat; raw_data + nidq_offsets(i)];
            end
            writematrix(nidq_concat, stitched_nidq_path);
            fprintf('Done.\n');
            
            % 4. Stitch IMEC
            fprintf('  > Stitching IMEC... ');
            imec_concat = [];
            for i = 1:nFiles
                fpath = fullfile(raw_imec(i).folder, raw_imec(i).name);
                raw_data = load(fpath);
                % Apply offset
                imec_concat = [imec_concat; raw_data + imec_offsets(i)];
            end
            writematrix(imec_concat, stitched_imec_path);
            fprintf('Done.\n');
            
            % Return new paths
            nidq_1hz_path = stitched_nidq_path;
            imec_1hz_path = stitched_imec_path;
        end
    end
end