
clc;
clear;

app = actxGetRunningServer('STK11.Application');
root = app.Personality2;

satNames = {'AutoSat1','AutoSat2','AutoSat3'};

for k = 1:numel(satNames)

    fprintf('\n===== %s =====\n', satNames{k});

    sat = root.GetObjectFromPath( ...
        ['*/Satellite/' satNames{k}]);

    % The COM interface identifies the propagator type
    prop = sat.Propagator;
    fprintf('Propagator interface: %s\n', class(prop));

    % Try reading the initial orbital elements
    try
        elems = prop.InitialState.Representation;

        fprintf('Initial state representation: %s\n', ...
            class(elems));

        try
            fprintf('Semi-major axis: %g\n', elems.SemiMajorAxis);
            fprintf('Eccentricity: %g\n', elems.Eccentricity);
            fprintf('Inclination: %g\n', elems.Inclination);
        catch
            fprintf(['Representation accessed; direct orbital ' ...
                     'element properties need further inspection.\n']);
        end

    catch ME
        fprintf('Orbital element query: %s\n', ME.message);
    end

end