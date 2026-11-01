%% ========================================================================
%% MODULE 08: VALIDATION OF MODULE 6 GEOMETRY AND MODULE 7 MEASUREMENTS
%% ========================================================================
% Needs: Module05_SatelliteStates.csv, Module06_EmitterGeometry.csv,
%        Module07_SyntheticMeasurements.csv
% Output: Module08_ValidationReport.csv
clc; format short g;

%% CONFIG (must match Modules 6 and 7)
cfg.m5file   = 'Module05_SatelliteStates.csv';
cfg.m6file   = 'Module06_EmitterGeometry.csv';
cfg.m7file   = 'Module07_SyntheticMeasurements.csv';
cfg.emitLat  = 31.00;  cfg.emitLon = 73.50;  cfg.emitAlt = 0.2;   % Module 6 truth
cfg.fc_Hz    = 2.4e9;      % Module 7 carrier
cfg.sigT_s   = 1e-6;       % Module 7 TDOA sigma
cfg.sigF_Hz  = 1.0;        % Module 7 FDOA sigma
cfg.fdoaSign = -1;         % Module 7: FDOA = -(fc/c)*dRR
cfg.refSat   = "AutoSat1";
cfg.minElDeg = 5;
cfg.altSat_km = 650;       % Module 2 orbit altitude
c    = 299792.458;         % km/s
mu   = 398600.4418;        % km^3/s^2
sats = ["AutoSat1","AutoSat2","AutoSat3"];
res  = struct('name',{},'pass',{},'info',{});

fprintf('\n===== MODULE 08 VALIDATION =====\n');

%% 1. Load Modules 5, 6 and 7
for f = {cfg.m5file, cfg.m6file, cfg.m7file}
    if ~isfile(f{1}), error('%s not found in current folder.', f{1}); end
end
S = readtable(cfg.m5file, setvartype(detectImportOptions(cfg.m5file, ...
    'VariableNamingRule','preserve'),'Epoch','string'));
G = readtable(cfg.m6file, setvartype(detectImportOptions(cfg.m6file, ...
    'VariableNamingRule','preserve'),'Epoch','string'));
o = detectImportOptions(cfg.m7file,'VariableNamingRule','preserve');
o = setvartype(o,{'Epoch','ReferenceSatellite','Satellite'},'string');
M = readtable(cfg.m7file, o);

S.Epoch = strtrim(string(S.Epoch));  S.Satellite = strtrim(string(S.Satellite));
G.Epoch = strtrim(string(G.Epoch));  G.Satellite = strtrim(string(G.Satellite));
M.Epoch = strtrim(M.Epoch);  M.Satellite = strtrim(M.Satellite);
M.ReferenceSatellite = strtrim(M.ReferenceSatellite);

[ue,~,~] = unique(G.Epoch,'stable');  ne = numel(ue);
fprintf('Loaded: %d M5 rows, %d M6 rows, %d M7 rows, %d epochs\n', ...
    height(S), height(G), height(M), ne);

%% 2. Independent recomputation of range, range rate and elevation
R = nan(ne,3);  RR = nan(ne,3);  EL = nan(ne,3);
maxDr = 0;  maxDv = 0;
for e = 1:ne
    for k = 1:3
        ig = find(G.Epoch==ue(e) & G.Satellite==sats(k), 1);
        is = find(S.Epoch==ue(e) & S.Satellite==sats(k), 1);
        if isempty(ig) || isempty(is)
            error('Missing row for %s at %s.', sats(k), ue(e));
        end
        rs = [S.X(is) S.Y(is) S.Z(is)];     vs = [S.VX(is) S.VY(is) S.VZ(is)];
        re = [G.Emit_X(ig) G.Emit_Y(ig) G.Emit_Z(ig)];
        ve = [G.Emit_VX(ig) G.Emit_VY(ig) G.Emit_VZ(ig)];
        dr = re - rs;   dv = ve - vs;
        R(e,k)  = norm(dr);
        RR(e,k) = dot(dr,dv)/norm(dr);
        EL(e,k) = asind(dot(-dr, re/norm(re))/norm(dr));  % emitter -> satellite vs local up  
        maxDr = max(maxDr, abs(R(e,k)  - G.Range_km(ig)));
        maxDv = max(maxDv, abs(RR(e,k) - G.RangeRate_kms(ig)));
    end
end
p = maxDr < 1e-6 && maxDv < 1e-9;
res(end+1) = report('Range / range-rate rebuilt from vectors', p, ...
    sprintf('max |dRange|=%.1e km, |dRR|=%.1e km/s', maxDr, maxDv));

%% 3. Earth-rotation velocity and emitter invariants
Ee = [G.Emit_X G.Emit_Y G.Emit_Z];   Ev = [G.Emit_VX G.Emit_VY G.Emit_VZ];
omegaE = [0 0 7.2921150e-5];
vExp = cross(repmat(omegaE, height(G), 1), Ee, 2);
eRot = max(abs(Ev(:) - vExp(:)));
a = 6378.137;  fl = 1/298.257223563;  e2 = fl*(2-fl);
phi = deg2rad(cfg.emitLat);  Np = a/sqrt(1 - e2*sin(phi)^2);
rxyExp = (Np + cfg.emitAlt)*cos(phi);
zExp   = (Np*(1-e2) + cfg.emitAlt)*sin(phi);
eInv = max([abs(hypot(Ee(:,1),Ee(:,2)) - rxyExp); abs(Ee(:,3) - zExp)]);
speedErr = max(abs(vecnorm(Ev,2,2) - 7.2921150e-5*rxyExp));
p = eRot < 1e-9 && eInv < 1e-6 && speedErr < 1e-9;
res(end+1) = report('Emitter v = omega x r; |xy| and z preserved', p, ...
    sprintf('dv=%.1e km/s, dpos=%.1e km, speed=%.4f km/s', ...
    eRot, eInv, mean(vecnorm(Ev,2,2))));

%% 4. Units and frame sanity for the satellites
rS = vecnorm([S.X S.Y S.Z],2,2);   vS = vecnorm([S.VX S.VY S.VZ],2,2);
rExp = a + cfg.altSat_km;
p = all(abs(rS - rExp) < 20) && max(abs(vS - sqrt(mu./rS))) < 0.05;
res(end+1) = report('Sat radius ~6378+650 km, speed = circular (km, km/s)', p, ...
    sprintf('r=%.1f..%.1f km, max|v-vcirc|=%.3f km/s', min(rS), max(rS), ...
    max(abs(vS - sqrt(mu./rS)))));

%% 5. Simultaneous visibility of the hidden emitter
p = all(EL(:) >= cfg.minElDeg);
res(end+1) = report(sprintf('All 3 sats above %d deg at every epoch', cfg.minElDeg), ...
    p, sprintf('min elevation=%.1f deg (per-sat min: %s, per-sat max: %s)', ...
    min(EL(:)), mat2str(min(EL,[],1),3), mat2str(max(EL,[],1),3)));
if ~p
    fprintf(['         Module 4 windows were computed for the VisibilityTestPoint ' ...
             '(31.52N, 74.36E), not the hidden emitter.\n']);
end
%% 6. Optional STK cross-check of range (frame consistency)
try
    app  = actxGetRunningServer('STK11.Application');
    scen = app.Personality2.CurrentScenario;
    em   = scen.Children.Item('HiddenEmitter');
    dR = nan(1,3);
    for k = 1:3
        ac = scen.Children.Item(char(sats(k))).GetAccessToObject(em);
        ac.ComputeAccess;
        pr = ac.DataProviders.Item('AER Data').Group.Item('Default');
        r  = pr.ExecSingle(char(ue(1)));
        d  = r.DataSets.GetDataSetByName('Range').GetValues;
        while iscell(d) && numel(d)==1, d = d{1}; end
        if iscell(d), d = d{1}; end
        dR(k) = str2double(string(d)) - R(1,k);
    end
    fprintf('  [INFO] STK minus MATLAB range, first epoch (km): %s\n', mat2str(dR,4));
    fprintf('         Tens of km is the expected GMST-only vs STK ICRF offset.\n');
catch ME
    fprintf('  [SKIP] STK range cross-check: %s\n', ME.message);
end

%% 7. Module 7: mapping, truth recompute, sign convention
[tf,loc] = ismember(M.Epoch, ue);
[tk,kk]  = ismember(M.Satellite, sats);
p = all(tf) && all(tk) && all(M.ReferenceSatellite == cfg.refSat) && ...
    height(M) == 2*ne && ~any(M.Satellite == cfg.refSat);
res(end+1) = report('M7 rows map to M6 epochs; ref=AutoSat1; 2 rows per epoch', p, ...
    sprintf('%d rows, %d epochs', height(M), ne));
if ~(all(tf) && all(tk)), error('Cannot continue: M7 epochs or satellites do not match M6.'); end

ik   = sub2ind(size(R), loc, kk);
dRng = R(ik)  - R(loc,1);
dRR  = RR(ik) - RR(loc,1);
T    = dRng/c;
F    = cfg.fdoaSign*cfg.fc_Hz/c*dRR;

eR = max(abs(dRng - M.RangeDifference_km));
eRRd = max(abs(dRR - M.RangeRateDifference_kms));
eT = max(abs(T - M.True_TDOA_s));
eF = max(abs(F - M.True_FDOA_Hz));
eFopp = max(abs(-F - M.True_FDOA_Hz));
p = eR < 1e-6 && eRRd < 1e-9 && eT < 1e-12 && eF < 1e-3;
res(end+1) = report('True TDOA/FDOA match independent recompute', p, ...
    sprintf('dRange=%.1e km, dRR=%.1e km/s, dTDOA=%.1e s, dFDOA=%.1e Hz', ...
    eR, eRRd, eT, eF));

p = eF < 1e-3 && eFopp > 100*max(eF,1e-3);
res(end+1) = report('Doppler sign: FDOA = -(fc/c)*dRR (receding => negative)', p, ...
    sprintf('error with this sign=%.1e Hz, with opposite sign=%.1e Hz', eF, eFopp));

%% 8. Noise: additivity and statistics
eSum = max(abs([M.Measured_TDOA_s - M.True_TDOA_s - M.Noise_TDOA_s; ...
                M.Measured_FDOA_Hz - M.True_FDOA_Hz - M.Noise_FDOA_Hz]));
p = eSum < 1e-9;
res(end+1) = report('Measured = True + Noise', p, sprintf('max residual=%.1e', eSum));

nT = M.Noise_TDOA_s;  nF = M.Noise_FDOA_Hz;  n = numel(nT);
tolMean = 4/sqrt(n);
pT = abs(mean(nT))/cfg.sigT_s < tolMean && abs(std(nT)/cfg.sigT_s - 1) < 0.5;
pF = abs(mean(nF))/cfg.sigF_Hz < tolMean && abs(std(nF)/cfg.sigF_Hz - 1) < 0.5;
res(end+1) = report('TDOA noise ~ N(0, sigT)', pT, ...
    sprintf('mean=%.2e s, std=%.2e s (cfg %.1e), N=%d', mean(nT), std(nT), cfg.sigT_s, n));
res(end+1) = report('FDOA noise ~ N(0, sigF)', pF, ...
    sprintf('mean=%.3f Hz, std=%.3f Hz (cfg %.1f), N=%d', mean(nF), std(nF), cfg.sigF_Hz, n));
fprintf('         Note: N=%d is small, so these statistics are only indicative.\n', n);

%% 9. Measurement quality context for Module 9
fprintf('  [INFO] |TDOA| range: %.2e to %.2e s (sigma %.1e s)\n', ...
    min(abs(M.True_TDOA_s)), max(abs(M.True_TDOA_s)), cfg.sigT_s);
fprintf('  [INFO] |FDOA| range: %.0f to %.0f Hz (sigma %.1f Hz)\n', ...
    min(abs(M.True_FDOA_Hz)), max(abs(M.True_FDOA_Hz)), cfg.sigF_Hz);
fprintf('  [INFO] 1-sigma TDOA = %.0f m range difference; 1-sigma FDOA = %.3f m/s\n', ...
    cfg.sigT_s*c*1000, cfg.sigF_Hz*c*1000/cfg.fc_Hz);

%% 10. Verdict and report
allPass = all([res.pass]);
fprintf('\n===== VERDICT =====\n%d/%d checks passed. ', sum([res.pass]), numel(res));
if allPass
    fprintf('Measurements are consistent: OK to proceed to Module 9.\n');
else
    fprintf('Fix the FAILED items above before Module 9.\n');
end
writetable(struct2table(res(:)), 'Module08_ValidationReport.csv');
fprintf('Report saved: Module08_ValidationReport.csv\nMODULE 08 FINISHED.\n');

%% Local functions (must stay at the end of the script)
function r = report(name, pass, info)
    if pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %s | %s\n', tag, name, info);
    r = struct('name', name, 'pass', pass, 'info', info);
end