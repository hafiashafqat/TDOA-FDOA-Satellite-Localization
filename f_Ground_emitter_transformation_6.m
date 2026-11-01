%% ========================================================================
%% MODULE 06: HIDDEN EMITTER + TRANSFORM TO SATELLITE FRAME
%% ========================================================================
% Needs: Satellite_State_Extraction.csv
% Output: Ground_emitter_transformation.csv
% Units: km, km/s. Frame: inertial (GMST-rotated) and satellite RIC.
clc; format short g;

%% 1. Connect to STK
try
    app = actxGetRunningServer('STK11.Application');
    root = app.Personality2;
    scenario = root.CurrentScenario;
catch ME
    error('Cannot connect to STK 11: %s', ME.message);
end
if isempty(scenario), error('No active STK scenario.'); end

%% 2. Create / reuse the hidden emitter (truth, not given to solver)
emitName = 'HiddenEmitter';
emitLat = 31.00;  emitLon = 73.50;  emitAlt = 0.2;   % deg, deg, km

if scenario.Children.Contains(8, emitName)
    emitter = scenario.Children.Item(emitName);
else
    emitter = scenario.Children.New(8, emitName);   % 8 = eFacility
end
emitter.Position.AssignGeodetic(emitLat, emitLon, emitAlt);
fprintf('Emitter "%s": lat %.3f, lon %.3f, alt %.3f km\n', ...
    emitName, emitLat, emitLon, emitAlt);

% WGS84 geodetic -> ECEF (km)
a = 6378.137; fl = 1/298.257223563; e2 = fl*(2-fl);
phi = deg2rad(emitLat); lam = deg2rad(emitLon);
N = a / sqrt(1 - e2*sin(phi)^2);
rECEF = [(N+emitAlt)*cos(phi)*cos(lam), ...
         (N+emitAlt)*cos(phi)*sin(lam), ...
         (N*(1-e2)+emitAlt)*sin(phi)];
omegaE = [0 0 7.2921150e-5];                 % Earth rate, rad/s

%% 3. Load Module 5 satellite states
f = 'Module05_SatelliteStates.csv';
if ~isfile(f), error('%s not found in current folder.', f); end
opts = detectImportOptions(f, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'Epoch', 'string');
S = readtable(f, opts);
S.Epoch = string(S.Epoch);

%% 4. Emitter inertial state at each unique epoch (ECEF -> inertial via GMST)
[ue, ~, idx] = unique(S.Epoch, 'stable');
E = zeros(numel(ue), 6);
for i = 1:numel(ue)
    t = datetime(char(ue(i)), 'InputFormat', 'dd MMM yyyy HH:mm:ss.SSS', ...
                 'Locale', 'en_US');
    jd   = datenum(t) + 1721058.5;                       % Julian date
    gmst = deg2rad(mod(280.46061837 + 360.98564736629*(jd - 2451545.0), 360));
    Rz = [cos(gmst) -sin(gmst) 0; sin(gmst) cos(gmst) 0; 0 0 1];
    re = (Rz*rECEF')';
    ve = cross(omegaE, re);                              % stationary on Earth
    E(i,:) = [re ve];
end

%% 5. Relative geometry and RIC transformation
n = height(S);
out = zeros(n, 20);
for r = 1:n
    rs = [S.X(r) S.Y(r) S.Z(r)];   vs = [S.VX(r) S.VY(r) S.VZ(r)];
    re = E(idx(r),1:3);            ve = E(idx(r),4:6);

    dr  = re - rs;                 % satellite -> emitter
    dv  = ve - vs;
    rho = norm(dr);
    rr  = dot(dr, dv) / rho;       % range rate, km/s

    Rh = rs / norm(rs);                        % radial
    Ch = cross(rs, vs); Ch = Ch / norm(Ch);    % cross-track
    Ih = cross(Ch, Rh);                        % in-track
    T  = [Rh; Ih; Ch];                         % inertial -> RIC

    out(r,:) = [re, ve, dr, dv, rho, rr, (T*dr')', (T*dv')'];
end

%% 6. Save
names = {'Emit_X','Emit_Y','Emit_Z','Emit_VX','Emit_VY','Emit_VZ', ...
         'dX','dY','dZ','dVX','dVY','dVZ','Range_km','RangeRate_kms', ...
         'dR_Radial','dR_InTrack','dR_CrossTrack', ...
         'dV_Radial','dV_InTrack','dV_CrossTrack'};
G = [S(:,{'Satellite','Epoch'}), array2table(out, 'VariableNames', names)];
writetable(G, 'Module06_EmitterGeometry.csv');

%% 7. Validation
fprintf('\n===== MODULE 06 VALIDATION =====\n');
fprintf('Records: %d | non-finite rows: %d\n', n, sum(any(~isfinite(out),2)));
fprintf('Emitter speed (expect ~0.399 km/s): %.4f km/s\n', norm(E(1,4:6)));
fprintf('Emitter radius (expect ~6371 km): %.2f km\n', norm(E(1,1:3)));
fprintf('Range min/max: %.1f / %.1f km\n', min(out(:,13)), max(out(:,13)));
disp(G(1:min(6,n), {'Satellite','Epoch','Range_km','RangeRate_kms'}));
fprintf('MODULE 06 FINISHED.\n');