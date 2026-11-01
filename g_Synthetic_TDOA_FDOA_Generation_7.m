
%% MODULE 7: TDOA / FDOA MEASUREMENT GENERATION
% Purpose:
% 1. Read satellite-emitter geometry from Module 6.
% 2. Calculate relative propagation-time differences (TDOA).
% 3. Calculate relative Doppler frequency differences (FDOA).
% 4. Add configurable Gaussian measurement noise.
% 5. Export synthetic measurements for Module 8.

clc;

fprintf('\n============================================\n');
fprintf(' MODULE 7: TDOA / FDOA MEASUREMENT GENERATION\n');
fprintf('============================================\n\n');

%% 1. CONFIGURATION PARAMETERS

inputFile  = 'Module06_EmitterGeometry.csv';
outputFile = 'Module07_SyntheticMeasurements.csv';

% RF configuration
c = 299792.458;              % Speed of light (km/s)
carrier_freq = 2.4e9;        % Assumed carrier frequency (Hz)

% Satellite configuration
referenceSatellite = "AutoSat1";

secondarySatellites = ["AutoSat2", "AutoSat3"];

% Measurement noise configuration
addNoise = true;

sigma_TDOA_s = 1e-6;         % TDOA noise standard deviation (seconds)
sigma_FDOA_Hz = 1;           % FDOA noise standard deviation (Hz)

% Reproducible random noise
rng(42, 'twister');

% Wavelength
lambda_km = c / carrier_freq;

fprintf('Carrier frequency: %.3f GHz\n', carrier_freq/1e9);
fprintf('Wavelength: %.6e km\n', lambda_km);
fprintf('Reference satellite: %s\n', referenceSatellite);

if addNoise
    fprintf('TDOA noise sigma: %.3f microseconds\n', ...
        sigma_TDOA_s*1e6);

    fprintf('FDOA noise sigma: %.3f Hz\n', sigma_FDOA_Hz);
else
    fprintf('Measurement noise: DISABLED\n');
end

%% 2. LOAD MODULE 6 OUTPUT

if ~isfile(inputFile)
    error('Module 6 output file not found: %s', inputFile);
end

geometry = readtable(inputFile);

requiredColumns = { ...
    'Satellite', ...
    'Epoch', ...
    'Range_km', ...
    'RangeRate_kms'};

missingColumns = setdiff(requiredColumns, ...
    geometry.Properties.VariableNames);

if ~isempty(missingColumns)
    error('Missing required columns: %s', ...
        strjoin(missingColumns, ', '));
end

% Normalize satellite and epoch identifiers
geometry.Satellite = strtrim(string(geometry.Satellite));
geometry.Epoch = strtrim(string(geometry.Epoch));

% Check numeric geometry
if ~isnumeric(geometry.Range_km) || ...
   ~isnumeric(geometry.RangeRate_kms)

    error('Range and range-rate columns must be numeric.');
end

if any(~isfinite(geometry.Range_km)) || ...
   any(~isfinite(geometry.RangeRate_kms))

    error('Module 6 contains invalid range or range-rate values.');
end

fprintf('\nModule 6 input loaded successfully.\n');
fprintf('Total geometry records: %d\n', height(geometry));

%% 3. IDENTIFY COMMON EPOCHS

epochs = unique(geometry.Epoch, 'stable');

if isempty(epochs)
    error('No valid epochs found in Module 6 output.');
end

satellitesRequired = [referenceSatellite, secondarySatellites];

nEpochs = numel(epochs);
nSecondary = numel(secondarySatellites);

fprintf('Unique measurement epochs: %d\n', nEpochs);

%% 4. PREALLOCATE OUTPUT ARRAYS

nMeasurements = nEpochs * nSecondary;

EpochOut = strings(nMeasurements, 1);
ReferenceOut = strings(nMeasurements, 1);
SatelliteOut = strings(nMeasurements, 1);

RangeDifference_km = zeros(nMeasurements, 1);
RangeRateDifference_kms = zeros(nMeasurements, 1);

True_TDOA_s = zeros(nMeasurements, 1);
True_FDOA_Hz = zeros(nMeasurements, 1);

Noise_TDOA_s = zeros(nMeasurements, 1);
Noise_FDOA_Hz = zeros(nMeasurements, 1);

Measured_TDOA_s = zeros(nMeasurements, 1);
Measured_FDOA_Hz = zeros(nMeasurements, 1);

%% 5. GENERATE TDOA AND FDOA MEASUREMENTS

row = 0;

for e = 1:nEpochs

    currentEpoch = epochs(e);

    % Extract records belonging to this epoch
    epochRows = geometry.Epoch == currentEpoch;

    epochData = geometry(epochRows, :);

    % Verify that each required satellite has exactly one record
    satelliteIndices = zeros(numel(satellitesRequired), 1);

    for s = 1:numel(satellitesRequired)

        matchingRows = find( ...
            epochData.Satellite == satellitesRequired(s));

        if numel(matchingRows) ~= 1
            error(['Expected exactly one record for %s at epoch %s. ' ...
                   'Found %d.'], ...
                   satellitesRequired(s), currentEpoch, ...
                   numel(matchingRows));
        end

        satelliteIndices(s) = matchingRows;
    end

    % Reference satellite geometry
    refIndex = satelliteIndices(1);

    R_ref = epochData.Range_km(refIndex);
    RR_ref = epochData.RangeRate_kms(refIndex);

    % Process each secondary satellite
    for s = 1:nSecondary

        row = row + 1;

        secondaryIndex = satelliteIndices(s + 1);

        R_secondary = epochData.Range_km(secondaryIndex);
        RR_secondary = epochData.RangeRate_kms(secondaryIndex);

        %% A. RANGE DIFFERENCE

        deltaR = R_secondary - R_ref;

        %% B. RANGE-RATE DIFFERENCE

        deltaRR = RR_secondary - RR_ref;

        %% C. TRUE TDOA

        trueTDOA = deltaR / c;

        %% D. TRUE FDOA

        trueFDOA = -(carrier_freq / c) * deltaRR;

        %% E. GENERATE CONFIGURABLE NOISE

        if addNoise

            noiseTDOA = sigma_TDOA_s * randn();

            noiseFDOA = sigma_FDOA_Hz * randn();

        else

            noiseTDOA = 0;
            noiseFDOA = 0;

        end

        %% F. NOISY MEASUREMENTS

        measuredTDOA = trueTDOA + noiseTDOA;

        measuredFDOA = trueFDOA + noiseFDOA;

        %% G. STORE RESULTS

        EpochOut(row) = currentEpoch;

        ReferenceOut(row) = referenceSatellite;

        SatelliteOut(row) = secondarySatellites(s);

        RangeDifference_km(row) = deltaR;

        RangeRateDifference_kms(row) = deltaRR;

        True_TDOA_s(row) = trueTDOA;

        True_FDOA_Hz(row) = trueFDOA;

        Noise_TDOA_s(row) = noiseTDOA;

        Noise_FDOA_Hz(row) = noiseFDOA;

        Measured_TDOA_s(row) = measuredTDOA;

        Measured_FDOA_Hz(row) = measuredFDOA;

    end

end

%% 6. CREATE OUTPUT TABLE

measurements = table( ...
    EpochOut, ...
    ReferenceOut, ...
    SatelliteOut, ...
    RangeDifference_km, ...
    RangeRateDifference_kms, ...
    True_TDOA_s, ...
    True_FDOA_Hz, ...
    Noise_TDOA_s, ...
    Noise_FDOA_Hz, ...
    Measured_TDOA_s, ...
    Measured_FDOA_Hz);

measurements.Properties.VariableNames = { ...
    'Epoch', ...
    'ReferenceSatellite', ...
    'Satellite', ...
    'RangeDifference_km', ...
    'RangeRateDifference_kms', ...
    'True_TDOA_s', ...
    'True_FDOA_Hz', ...
    'Noise_TDOA_s', ...
    'Noise_FDOA_Hz', ...
    'Measured_TDOA_s', ...
    'Measured_FDOA_Hz'};

%% 7. VALIDATE GENERATED MEASUREMENTS

numericData = measurements{:, 4:end};

if any(~isfinite(numericData), 'all')
    error('Invalid numerical values detected in Module 7 output.');
end

expectedRows = nEpochs * nSecondary;

if height(measurements) ~= expectedRows
    error('Unexpected number of measurement records.');
end

fprintf('\nMeasurement validation successful.\n');
fprintf('Expected measurement records: %d\n', expectedRows);
fprintf('Generated measurement records: %d\n', height(measurements));

%% 8. EXPORT RESULTS

writetable(measurements, outputFile);

fprintf('\nModule 7 completed successfully.\n');
fprintf('Output file: %s\n', fullfile(pwd, outputFile));

%% 9. DISPLAY SAMPLE RESULTS

fprintf('\nFirst measurement records:\n');

disp(measurements(1:min(10, height(measurements)), :));

fprintf('\n============================================\n');
fprintf(' MODULE 7 EXECUTION FINISHED\n');
fprintf('============================================\n');
