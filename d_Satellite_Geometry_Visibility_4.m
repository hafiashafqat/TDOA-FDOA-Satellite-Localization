%% ========================================================================
%% MODULE 04: SATELLITE GEOMETRY, ACCESS & COMMON VISIBILITY
%% ========================================================================
clc;
clearvars -except app root scenario; % Retain persistent COM connections if they exist
format short g;

%% 1. Connect to the existing STK session
try
    app = actxGetRunningServer('STK11.Application');
    root = app.Personality2;
    scenario = root.CurrentScenario;
catch ME
    error('Cannot connect to the existing STK 11 session: %s', ME.message);
end

if isempty(scenario)
    error('No scenario is currently open in STK.');
end

fprintf('\nExecuting Module 4: Access & Common Visibility Integration\n');
fprintf('Active Scenario: %s\n', scenario.InstanceName);

%% 2. Check and Validate Active Satellites
satNames = {'AutoSat1','AutoSat2','AutoSat3'};
satObjects = cell(1,3);

for k = 1:3
    try
        satObjects{k} = scenario.Children.Item(satNames{k});
        fprintf('  Satellite validated: %s\n', satNames{k});
    catch
        error('Critical Failure: Satellite %s is missing from the tree.', satNames{k});
    end
end

%% 3. Target facility = the hidden emitter (position always enforced)
facilityName = 'HiddenEmitter';
targetLat = 31.00;
targetLon = 73.50;
targetAlt = 0.2;      % km, same as Module 6

if scenario.Children.Contains(8, facilityName)
    facility = scenario.Children.Item(facilityName);
    fprintf('  Facility "%s" detected in scenario tree.\n', facilityName);
else
    fprintf('  Facility "%s" missing. Creating it...\n', facilityName);
    facility = scenario.Children.New(8, facilityName);   % 8 = eFacility
end
facility.Position.AssignGeodetic(targetLat, targetLon, targetAlt);
fprintf('  Facility at lat %.3f, lon %.3f, alt %.3f km\n', ...
    targetLat, targetLon, targetAlt);
%% 4. Verify/Enforce the 5-Degree Minimum Elevation Constraint
try
    % Access the base constraints collection
    accessConstraints = facility.AccessConstraints;
    
    % FIX: Check if constraint is already active, otherwise add it safely
    try
        constraint = accessConstraints.AddConstraint(14); % 14 = eCstrElevation
    catch ME
        if contains(ME.message, 'Constraint already active')
            % Retrieve the existing constraint handle instead
            constraint = accessConstraints.GetActiveConstraint(14);
        else
            rethrow(ME);
        end
    end

    % Explicitly activate and enforce the constraint properties
    constraint.EnableMin = 1; 
    constraint.Min = 5.0; 

    fprintf('  Constraint Enforced: %s Minimum Elevation = %.1f deg\n', ...
        facilityName, double(constraint.Min));
catch ME
    error('Cannot verify or enforce the elevation constraint. Details: %s', ME.message);
end


%% 5. Extract Satellite Positions at Scenario Reference Epoch
positionRows = cell(0,5);
startTimeStr = scenario.StartTime;

for k = 1:3
    try
        % Extract positional states in ICRF system frame
        provider = satObjects{k}.DataProviders.Item('Cartesian Position').Group.Item('ICRF');
        result = provider.ExecSingle(startTimeStr);
        
        x = scalarValue(result.DataSets.GetDataSetByName('x').GetValues);
        y = scalarValue(result.DataSets.GetDataSetByName('y').GetValues);
        z = scalarValue(result.DataSets.GetDataSetByName('z').GetValues);

        positionRows(end+1,:) = {satNames{k}, char(startTimeStr), x, y, z};
        fprintf('    %s Base Position -> X: %g, Y: %g, Z: %g km\n', satNames{k}, x, y, z);
    catch ME
        error('Position extraction failed for %s: %s', satNames{k}, ME.message);
    end
end

positionTable = cell2table(positionRows, ...
    'VariableNames', {'Satellite','Epoch','X_km','Y_km','Z_km'});
writetable(positionTable, 'Module04_SatellitePositions.csv');

%% 6. Compute Individual Access Intervals
accessTables = cell(1,3);

for k = 1:3
    fprintf('  Computing analytical access for %s...\n', satNames{k});
    try
        access = satObjects{k}.GetAccessToObject(facility);
        access.ComputeAccess();

        provider = access.DataProviders.Item('Access Data');
        result = provider.Exec(scenario.StartTime, scenario.StopTime);

        starts = normalizeValues(result.DataSets.GetDataSetByName('Start Time').GetValues);
        stops = normalizeValues(result.DataSets.GetDataSetByName('Stop Time').GetValues);

        if isempty(starts)
            starts = {};
            stops = {};
        end

        accessTables{k} = table(...
            repmat(string(satNames{k}), numel(starts), 1), ...
            string(starts(:)), string(stops(:)), ...
            'VariableNames', {'Satellite','Start','Stop'});

        fprintf('    Intervals discovered: %d\n', numel(starts));
    catch ME
        error('Access loop failed for satellite %s: %s', satNames{k}, ME.message);
    end
end

allAccess = vertcat(accessTables{:});
writetable(allAccess, 'Module04_AccessIntervals.csv');

%% 7. Compute Overlapping Common Visibility (3-Way Sat Intersect)
common = accessTables{1}(:, {'Start','Stop'});

for k = 2:3
    next = accessTables{k};
    overlaps = cell(0,2);

    for i = 1:height(common)
        a = parseTime(common.Start(i));
        b = parseTime(common.Stop(i));

        for j = 1:height(next)
            c = parseTime(next.Start(j));
            d = parseTime(next.Stop(j));

            lo = max(a,c);
            hi = min(b,d);

            if hi > lo
                overlaps(end+1,:) = {string(formatTime(lo)), string(formatTime(hi))};
            end
        end
    end

    common = cell2table(overlaps, 'VariableNames', {'Start','Stop'});
end

%% 8. Export Metrics and Summarize
writetable(common, 'Module04_CommonVisibility.csv');

fprintf('\n========== MODULE 04 VERIFICATION SUMMARY ==========\n');
fprintf('  Target Facility Tracked    : %s\n', facilityName);
fprintf('  Position Vectors Processed : %d\n', height(positionTable));
fprintf('  Total Indiv. Access Records: %d\n', height(allAccess));
fprintf('  Intersection Multi-Access  : %d windows\n', height(common));
fprintf('===================================================\n');
disp(common);

%% Local helper functions
function value = scalarValue(data)
    while iscell(data) && numel(data) == 1
        data = data{1};
    end
    value = double(data);
    value = value(1);
end

function values = normalizeValues(data)
    if isempty(data), values = {}; return; end
    while iscell(data) && numel(data) == 1 && iscell(data{1})
        data = data{1};
    end
    if ischar(data) || isstring(data)
        data = cellstr(data);
    elseif ~iscell(data)
        data = num2cell(data);
    end
    values = cellfun(@(v) char(string(v)), data(:), 'UniformOutput', false);
end

function t = parseTime(value)
    formats = {'d MMM yyyy HH:mm:ss.SSS', 'd MMM yyyy HH:mm:ss'};
    for k = 1:numel(formats)
        try
            t = datetime(char(value), 'InputFormat', formats{k}, 'Locale', 'en_US');
            return;
        catch
        end
    end
    error('Could not parse STK time string: %s', char(value));
end
function s = formatTime(t)
    % Set the display format property directly on the datetime object
    t.Format = 'd MMM yyyy HH:mm:ss.SSS';
    s = char(t);
end

