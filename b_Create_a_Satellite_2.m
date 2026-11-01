%% ========================================================================
%% MODULE 2: CREATE AND PROPAGATE SATELLITES
%% ========================================================================
clc;

% Connect to the existing STK session and scenario
app = actxGetRunningServer('STK11.Application');
root = app.Personality2;
scenario = root.CurrentScenario;

% Definitive constellation parameters
names = {'AutoSat1','AutoSat2','AutoSat3'};
phase = [0, 10, 20];         % True anomaly (degrees) along track
raan_offsets = [0, 5, 10];   % Spreads planes out so they aren't collinear
altitude = 650;              % Circular orbit altitude (km)
inclination = 60;            % Inclination (degrees)

for k = 1:numel(names)
    
    % Defensive check: If satellite exists from a previous run, delete it cleanly
    if scenario.Children.Contains(18, names{k})
        fprintf('Satellite %s already exists. Resetting object...\n', names{k});
        existingSat = scenario.Children.Item(names{k});
        existingSat.Unload();
    end

    % Create unique satellite instance (18 = eSatellite object type)
    sat = scenario.Children.New(18, names{k});

    % Select analytical two-body orbit propagation
    sat.SetPropagatorType('ePropagatorTwoBody');
    prop = sat.Propagator;

    % Configure classical orbital elements interface
    oe = prop.InitialState.Representation.ConvertTo('eOrbitStateClassical');

    oe.SizeShapeType = 'eSizeShapeAltitude';
    oe.SizeShape.PerigeeAltitude = altitude;
    oe.SizeShape.ApogeeAltitude = altitude;

    oe.Orientation.Inclination = inclination;
    oe.Orientation.AscNode.Value = raan_offsets(k); % Spatially separates orbital planes

    oe.LocationType = 'eLocationTrueAnomaly';
    oe.Location.Value = phase(k);

    % Commit changes to the STK Engine and propagate the path
    prop.InitialState.Representation.Assign(oe);
    prop.Propagate;

    fprintf(' %s initialized and successfully propagated.\n', names{k});
end

disp('MODULE 02 COMPLETE');

