
%% ========================================================================
%% MODULE 05: SATELLITE STATE EXTRACTION
%% ========================================================================
% STK 11 + MATLAB
% Uses Module04_CommonVisibility.csv
% Extracts satellite position and velocity in the ICRF provider group.
%
% Output:
%   Module05_SatelliteStates.csv
%
% Important:
%   - Reuses the existing STK scenario and satellites.
%   - Does not create, delete, or modify STK objects.
%   - Samples common visibility intervals at 60-second intervals.
%   - Includes each interval's start time.
%   - Position units: km
%   - Velocity units: km/s

clc;
format short g;

fprintf('\n====================================================\n');
fprintf(' MODULE 05: SATELLITE STATE EXTRACTION\n');
fprintf('====================================================\n');

%% 1. Connect to the existing STK session
try
    app = actxGetRunningServer('STK11.Application');
    root = app.Personality2;
    scenario = root.CurrentScenario;
catch ME
    error('Cannot connect to the running STK 11 session: %s', ...
        ME.message);
end

if isempty(scenario)
    error('No active STK scenario is open.');
end

fprintf('Active scenario: %s\n', scenario.InstanceName);

%% 2. Verify the three existing satellites
satNames = {'AutoSat1','AutoSat2','AutoSat3'};
satObjects = cell(1,3);

for k = 1:3
    try
        satObjects{k} = scenario.Children.Item(satNames{k});
        fprintf('Satellite verified: %s\n', satNames{k});
    catch
        error('Required satellite "%s" was not found.', satNames{k});
    end
end


%% 3. Load Module 4 common visibility intervals

visibilityFile = fullfile(pwd, 'Module04_CommonVisibility.csv');

if ~isfile(visibilityFile)
    error('Visibility file not found in current folder: %s', visibilityFile);
end

fprintf('Reading visibility file: %s\n', visibilityFile);

% First try the normal table reader.
visibility = readtable(visibilityFile, 'VariableNamingRule', 'preserve');

% Repair column names if MATLAB imported them differently.
names = string(visibility.Properties.VariableNames);
namesLower = lower(strtrim(names));

startCol = find(namesLower == "start", 1);
stopCol  = find(namesLower == "stop", 1);

if isempty(startCol) || isempty(stopCol)
    % Fall back to reading the first two CSV fields as text.
    raw = readcell(visibilityFile);

    if size(raw, 2) < 2 || size(raw, 1) < 2
        error('Visibility CSV does not contain usable interval data.');
    end

    % Identify Start/Stop columns from header names, if present.
    header = lower(strtrim(string(raw(1,:))));
    startCol = find(header == "start", 1);
    stopCol  = find(header == "stop", 1);

    if isempty(startCol) || isempty(stopCol)
        % The saved file may have been split into multiple date fields.
        % Stop rather than silently interpreting the wrong columns.
        error(['Could not identify Start and Stop columns in %s. ' ...
               'Please send the output of: disp(readcell(''%s''))'], ...
               visibilityFile, visibilityFile);
    end

    visibility = table( ...
        string(raw(2:end,startCol)), ...
        string(raw(2:end,stopCol)), ...
        'VariableNames', {'Start','Stop'});
else
    visibility = table( ...
        string(visibility{:,startCol}), ...
        string(visibility{:,stopCol}), ...
        'VariableNames', {'Start','Stop'});
end

% Remove blank rows and validate timestamps.
keep = strlength(strtrim(visibility.Start)) > 0 & ...
       strlength(strtrim(visibility.Stop)) > 0;
visibility = visibility(keep,:);

if isempty(visibility)
    error('No valid visibility intervals were found.');
end

fprintf('Visibility intervals loaded: %d\n', height(visibility));
disp(visibility);


%% 4. Generate observation epochs at 60-second cadence
% Windows are trimmed at both ends so epochs stay clear of the 5 deg edge
% (margin for the STK vs GMST-only emitter offset, about 0.5 deg).

cadenceSeconds = 60;
trimSeconds    = 30;
epochs = datetime.empty(0,1);

for i = 1:height(visibility)
    intervalStart = parseSTKTime(visibility.Start(i)) + seconds(trimSeconds);
    intervalStop  = parseSTKTime(visibility.Stop(i))  - seconds(trimSeconds);

    if intervalStop <= intervalStart
        warning('Interval %d is too short after trimming; skipping.', i);
        continue;
    end

    intervalEpochs = (intervalStart:seconds(cadenceSeconds):intervalStop)';
    epochs = [epochs; intervalEpochs]; %#ok<AGROW>
end

epochs = sort(unique(epochs));

if isempty(epochs)
    error('No valid observation epochs were generated.');
end

fprintf('Trim per edge: %d s | Sampling cadence: %d s\n', trimSeconds, cadenceSeconds);
fprintf('Total observation epochs: %d\n', numel(epochs));
fprintf('Expected satellite-state records: %d\n', ...
    numel(epochs) * numel(satNames));

%% 5. Prepare ICRF position and velocity providers
positionProviders = cell(1,3);
velocityProviders = cell(1,3);

for k = 1:3
    try
        positionProviders{k} = ...
            satObjects{k}.DataProviders.Item( ...
            'Cartesian Position').Group.Item('ICRF');

        velocityProviders{k} = ...
            satObjects{k}.DataProviders.Item( ...
            'Cartesian Velocity').Group.Item('ICRF');

        fprintf('ICRF providers ready: %s\n', satNames{k});
    catch ME
        error(['Could not access the ICRF position/velocity providers ' ...
            'for %s: %s'], satNames{k}, ME.message);
    end
end

%% 6. Extract position and velocity at every epoch
% Each row contains one satellite at one observation epoch.

nExpected = numel(epochs) * numel(satNames);
stateRows = cell(nExpected, 8);
row = 0;

for e = 1:numel(epochs)
    epochString = datestr(epochs(e), 'dd mmm yyyy HH:MM:SS.FFF');

    fprintf('Epoch %d/%d: %s\n', ...
        e, numel(epochs), epochString);

    for k = 1:numel(satNames)
        try
            posResult = positionProviders{k}.ExecSingle(epochString);
            velResult = velocityProviders{k}.ExecSingle(epochString);

            x  = getScalar(posResult, 'x');
            y  = getScalar(posResult, 'y');
            z  = getScalar(posResult, 'z');

            vx = getScalar(velResult, 'x');
            vy = getScalar(velResult, 'y');
            vz = getScalar(velResult, 'z');

            values = [x y z vx vy vz];

            if any(~isfinite(values))
                error('Non-finite position or velocity value returned.');
            end

            row = row + 1;
            stateRows(row,:) = { ...
                satNames{k}, epochString, x, y, z, vx, vy, vz};

        catch ME
            error(['State extraction failed for %s at %s: %s'], ...
                satNames{k}, epochString, ME.message);
        end
    end
end

%% 7. Create and save the output table
stateRows = stateRows(1:row,:);

states = cell2table(stateRows, ...
    'VariableNames', ...
    {'Satellite','Epoch','X','Y','Z','VX','VY','VZ'});

outputFile = 'Module05_SatelliteStates.csv';
writetable(states, outputFile);

fprintf('\nSaved satellite states: %s\n', outputFile);

%% 8. Validate extracted records
fprintf('\n========== MODULE 05 VALIDATION ==========\n');

fprintf('Total records: %d\n', height(states));
fprintf('Unique epochs: %d\n', numel(unique(states.Epoch)));
fprintf('Expected records: %d\n', ...
    numel(epochs) * numel(satNames));

numericData = states{:, {'X','Y','Z','VX','VY','VZ'}};
missingCount = sum(any(~isfinite(numericData), 2));

fprintf('Rows containing non-finite values: %d\n', missingCount);

[epochGroups, ~, groupIndex] = unique(states.Epoch);
counts = accumarray(groupIndex, 1);
badEpochCount = sum(counts ~= numel(satNames));

fprintf('Epochs with exactly 3 satellite records: %d/%d\n', ...
    sum(counts == numel(satNames)), numel(epochGroups));

fprintf('Epochs with an incorrect satellite count: %d\n', ...
    badEpochCount);

if height(states) ~= nExpected || missingCount > 0 || badEpochCount > 0
    warning('Module 05 output did not pass every validation check.');
else
    fprintf('Record-count and numeric validation passed.\n');
end

fprintf('Position units: km\n');
fprintf('Velocity units: km/s\n');
fprintf('Provider group: ICRF\n');

disp(states(1:min(9,height(states)),:));

fprintf('\nMODULE 05 FINISHED.\n');

%% ========================================================================
%% LOCAL HELPER FUNCTIONS
%% ========================================================================

function t = parseSTKTime(value)
    % Parse timestamps written by Module 4.
    formats = { ...
        'dd MMM yyyy HH:mm:ss.SSS', ...
        'd MMM yyyy HH:mm:ss.SSS', ...
        'dd MMM yyyy HH:mm:ss', ...
        'd MMM yyyy HH:mm:ss'};

    for n = 1:numel(formats)
        try
            t = datetime(char(value), ...
                'InputFormat', formats{n}, ...
                'Locale', 'en_US');
            return;
        catch
        end
    end

    error('Unable to parse visibility timestamp: %s', char(value));
end

function value = getScalar(result, dataSetName)
    % Read a scalar value from an STK data-provider result.
    data = result.DataSets.GetDataSetByName(dataSetName).GetValues;

    % Unwrap nested cell arrays returned by some STK COM interfaces.
    while iscell(data) && numel(data) == 1
        data = data{1};
    end

    if iscell(data)
        data = data{1};
    end

    value = str2double(string(data));

    if ~isscalar(value) || ~isfinite(value)
        error('Invalid scalar returned for dataset "%s".', dataSetName);
    end
end
