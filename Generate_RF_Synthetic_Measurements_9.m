
%% Generate_RF_Synthetic_Measurements.m
% Purpose:
% Generate synthetic TDOA and FDOA observations for three STK satellites.
% Satellite states are taken from RF_Measurement_Window.mat.
% The Facility Cartesian Position COM provider is NOT used.
%
% Outputs:
%   RF_Synthetic_Measurements.mat  - measurements for the estimator
%   RF_Validation_Truth.mat        - truth coordinates for validation only

clearvars;
clc;

fprintf('\n============================================\n');
fprintf(' SYNTHETIC RF MEASUREMENT GENERATION\n');
fprintf('============================================\n');

%% 1. File checks
stateFile  = 'RF_Measurement_Window.mat';
configFile = 'RF_Localization_Config.mat';

if ~isfile(stateFile)
    error('Cannot find %s in the current folder.', stateFile);
end

if ~isfile(configFile)
    error('Cannot find %s in the current folder.', configFile);
end

S = load(stateFile);
C = load(configFile);

if ~isfield(S, 'measurementWindow')
    error(['RF_Measurement_Window.mat does not contain ', ...
           'measurementWindow. Inspect its actual variable names.']);
end

if ~isfield(C, 'CFG')
    error('RF_Localization_Config.mat does not contain CFG.');
end

measurementWindow = S.measurementWindow;
CFG = C.CFG;

%% 2. Measurement settings
numberOfSatellites = 3;
sampleStep_s = 10;
f0_Hz = CFG.mission.targetFrequency_Hz;
sigmaTDOA_s = CFG.measurement.TDOA_noiseStd_s;
sigmaFDOA_Hz = CFG.measurement.FDOA_noiseStd_Hz;

% Observation window already validated in the STK scenario.
% There are 33 samples spanning 320 seconds.
startTime = datetime('06 Oct 2026 20:12:47.000', ...
    'InputFormat', 'dd MMM yyyy HH:mm:ss.SSS', ...
    'TimeZone', 'UTC');

%% 3. Load satellite positions and velocities
if ~isfield(measurementWindow, 'state')
    error('measurementWindow does not contain the expected state field.');
end

if numel(measurementWindow.state) ~= numberOfSatellites
    error('Expected exactly three satellite state records.');
end

rSat_km = cell(numberOfSatellites, 1);
vSat_kmps = cell(numberOfSatellites, 1);

for k = 1:numberOfSatellites
    state = measurementWindow.state(k);

    if ~isfield(state, 'position_km') || ...
       ~isfield(state, 'velocity_kmps')
        error(['Satellite state %d does not contain position_km ', ...
               'and velocity_kmps.'], k);
    end

    rSat_km{k} = toThreeByN(state.position_km);
    vSat_kmps{k} = toThreeByN(state.velocity_kmps);

    if ~isequal(size(rSat_km{k}), size(vSat_kmps{k}))
        error('Position and velocity dimensions differ for satellite %d.', k);
    end
end

N = size(rSat_km{k}, 2);

for k = 1:numberOfSatellites
    if size(rSat_km{k}, 2) ~= N
        error('Satellite sample counts do not match.');
    end
end

if N ~= 33
    error('Expected 33 samples, but found %d. Check the observation window.', N);
end

if any(~isfinite(rSat_km{1}(:))) || ...
   any(~isfinite(rSat_km{2}(:))) || ...
   any(~isfinite(rSat_km{3}(:))) || ...
   any(~isfinite(vSat_kmps{1}(:))) || ...
   any(~isfinite(vSat_kmps{2}(:))) || ...
   any(~isfinite(vSat_kmps{3}(:)))
    error('Satellite state data contains NaN or Inf.');
end

fprintf('Satellite samples: %d per satellite\n', N);
fprintf('Sampling interval: %.1f seconds\n', sampleStep_s);
fprintf('Signal frequency: %.3f MHz\n', f0_Hz / 1e6);
fprintf('TDOA noise std: %.1f ns\n', sigmaTDOA_s * 1e9);
fprintf('FDOA noise std: %.2f Hz\n', sigmaFDOA_Hz);

%% 4. Read hidden validation location
truth = CFG.transmitter.validationTruth;

lat_deg = truth.latitude_deg;
lon_deg = truth.longitude_deg;
alt_m = truth.altitude_m;

% The truth coordinates are used only to synthesize observations.
% They are not included in RF_Synthetic_Measurements.mat.

%% 5. Convert ground location from WGS-84 to ECEF
% WGS-84 constants
a_km = 6378.137;
flattening = 1 / 298.257223563;
e2 = flattening * (2 - flattening);

lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);
h_km = alt_m / 1000;

Nphi = a_km / sqrt(1 - e2 * sin(lat)^2);

rECEF_km = [
    (Nphi + h_km) * cos(lat) * cos(lon);
    (Nphi + h_km) * cos(lat) * sin(lon);
    (Nphi * (1 - e2) + h_km) * sin(lat)
];

%% 6. Calculate ground emitter inertial state at each epoch
% Rotate ECEF coordinates into an Earth-centred inertial frame using
% Greenwich Mean Sidereal Time. This is an approximate frame conversion;
% it neglects polar motion and other small Earth-orientation corrections.
%
% This calculation is appropriate only if the satellite states use a
% compatible Earth-centred inertial frame and the timestamps are aligned.

c_kmps = 299792.458;
omegaEarth_radps = 7.2921150e-5;

times = startTime + seconds((0:N-1)' * sampleStep_s);

rGround_km = zeros(3, N);
vGround_kmps = zeros(3, N);

for n = 1:N
    jd = juliandate(times(n));
    T = (jd - 2451545.0) / 36525;

    gmst_deg = 280.46061837 ...
        + 360.98564736629 * (jd - 2451545.0) ...
        + 0.000387933 * T^2 ...
        - T^3 / 38710000;

    theta = deg2rad(mod(gmst_deg, 360));

    Rz = [
        cos(theta), -sin(theta), 0;
        sin(theta),  cos(theta), 0;
        0,           0,          1
    ];

    rGround_km(:, n) = Rz * rECEF_km;

    % Inertial velocity of a ground-fixed point due to Earth's rotation
    vGround_kmps(:, n) = cross( ...
        [0; 0; omegaEarth_radps], rGround_km(:, n));
end

%% 7. Calculate range and range rate for each satellite
range_km = zeros(numberOfSatellites, N);
rangeRate_kmps = zeros(numberOfSatellites, N);

for k = 1:numberOfSatellites
    for n = 1:N
        relativePosition = rSat_km{k}(:, n) - rGround_km(:, n);
        relativeVelocity = vSat_kmps{k}(:, n) - vGround_kmps(:, n);

        range_km(k, n) = norm(relativePosition);

        if range_km(k, n) <= 0
            error('Invalid range for satellite %d at sample %d.', k, n);
        end

        rangeRate_kmps(k, n) = ...
            dot(relativePosition, relativeVelocity) / range_km(k, n);
    end
end

%% 8. Generate noiseless TDOA and FDOA
% Satellite 1 is the reference.
%
% TDOA(i,1) = (range_i - range_1)/c
%
% Received frequency:
% f_received = f0 * (1 - rangeRate/c)
%
% FDOA(i,1) = receivedFrequency_i - receivedFrequency_1

TDOA_clean_s = zeros(N, numberOfSatellites - 1);
FDOA_clean_Hz = zeros(N, numberOfSatellites - 1);

receivedFrequency_Hz = f0_Hz .* ...
    (1 - rangeRate_kmps ./ c_kmps);

for k = 2:numberOfSatellites
    col = k - 1;

    TDOA_clean_s(:, col) = ...
        (range_km(k, :) - range_km(1, :))' / c_kmps;

    FDOA_clean_Hz(:, col) = ...
        (receivedFrequency_Hz(k, :) - ...
         receivedFrequency_Hz(1, :))';
end

%% 9. Add measurement noise
rng(42, 'twister');  % Reproducible synthetic test

TDOA_noise_s = sigmaTDOA_s * randn(N, 2);
FDOA_noise_Hz = sigmaFDOA_Hz * randn(N, 2);

TDOA_measured_s = TDOA_clean_s + TDOA_noise_s;
FDOA_measured_Hz = FDOA_clean_Hz + FDOA_noise_Hz;

%% 10. Package estimator input
% Columns correspond to:
%   column 1: Satellite 2 relative to Satellite 1
%   column 2: Satellite 3 relative to Satellite 1
%
% No emitter latitude/longitude is saved in this measurement structure.

measurements.time_UTC = times;
measurements.sampleStep_s = sampleStep_s;
measurements.referenceSatellite = 'AutoSat1';
measurements.otherSatellites = {'AutoSat2', 'AutoSat3'};
measurements.frequency_Hz = f0_Hz;

measurements.TDOA_s = TDOA_measured_s;
measurements.FDOA_Hz = FDOA_measured_Hz;

measurements.TDOA_noiseStd_s = sigmaTDOA_s;
measurements.FDOA_noiseStd_Hz = sigmaFDOA_Hz;

measurements.numberOfSamples = N;
measurements.measurementType = 'Synthetic';
measurements.frameAssumption = ...
    'Satellite states assumed compatible with the approximate inertial ground state';

% Save only the measurement structure for the future estimator.
save('RF_Synthetic_Measurements.mat', 'measurements');

%% 11. Save truth separately for validation only
validationTruth.latitude_deg = lat_deg;
validationTruth.longitude_deg = lon_deg;
validationTruth.altitude_m = alt_m;

save('RF_Validation_Truth.mat', 'validationTruth');

%% 12. Report results
fprintf('\n--------------------------------------------\n');
fprintf('MEASUREMENT GENERATION COMPLETE\n');
fprintf('--------------------------------------------\n');
fprintf('Number of epochs: %d\n', N);
fprintf('TDOA measurements generated: %d\n', numel(TDOA_measured_s));
fprintf('FDOA measurements generated: %d\n', numel(FDOA_measured_Hz));

fprintf('\nTDOA range: %.6f to %.6f microseconds\n', ...
    min(TDOA_measured_s(:)) * 1e6, ...
    max(TDOA_measured_s(:)) * 1e6);

fprintf('FDOA range: %.3f to %.3f Hz\n', ...
    min(FDOA_measured_Hz(:)), ...
    max(FDOA_measured_Hz(:)));

fprintf('\nCreated:\n');
fprintf('  RF_Synthetic_Measurements.mat\n');
fprintf('  RF_Validation_Truth.mat\n');

fprintf('\nNext step: validate these measurements before localization.\n');

%% Local helper: convert stored state data to a 3-by-N numeric matrix
function A = toThreeByN(raw)

    if iscell(raw)
        try
            A = str2double(string(raw));
        catch
            error('Could not convert stored state data from cell format.');
        end
    else
        A = double(raw);
    end

    if ~isnumeric(A) || isempty(A) || ~ismatrix(A)
        error('State data must be a nonempty numeric matrix.');
    end

    if size(A, 2) == 3
        A = A.';
    elseif size(A, 1) ~= 3
        error('State matrix must have three coordinates per sample.');
    end

    if any(~isfinite(A(:)))
        error('State matrix contains nonfinite values.');
    end
end