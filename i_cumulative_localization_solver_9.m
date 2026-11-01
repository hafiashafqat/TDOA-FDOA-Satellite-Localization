%% ========================================================================
%% MODULE 09: CUMULATIVE LOCALIZATION SOLVER (TDOA + FDOA)
%% ========================================================================
% Needs: Module05_SatelliteStates.csv, Module07_SyntheticMeasurements.csv
% Output: Module09_LocalizationResults.csv
clc; format short g;

%% CONFIG (must match Modules 6 and 7)
cfg.m5file = 'Module05_SatelliteStates.csv';
cfg.m7file = 'Module07_SyntheticMeasurements.csv';
cfg.fc_Hz = 2.4e9;  cfg.sigT_s = 1e-6;  cfg.sigF_Hz = 1.0;
cfg.fdoaSign = -1;
cfg.assumedAlt_km = 0.2;                 % emitter altitude assumed known
cfg.omegaE = 7.2921150e-5;
cfg.truth = [31.00 73.50 0.2];           % scoring ONLY, not used by solver
cfg.gridHalf = 12;  cfg.gridStep = 1;    % deg, search grid around sub-points
c = 299792.458;
sats = ["AutoSat1","AutoSat2","AutoSat3"];
cfg.c = c;

%% 1. Load data
S = readtable(cfg.m5file, setvartype(detectImportOptions(cfg.m5file, ...
    'VariableNamingRule','preserve'),'Epoch','string'));
o = detectImportOptions(cfg.m7file,'VariableNamingRule','preserve');
o = setvartype(o,{'Epoch','ReferenceSatellite','Satellite'},'string');
M = readtable(cfg.m7file, o);
S.Epoch = strtrim(string(S.Epoch));  S.Satellite = strtrim(string(S.Satellite));
M.Epoch = strtrim(M.Epoch);          M.Satellite = strtrim(M.Satellite);

ue = unique(S.Epoch,'stable');  ne = numel(ue);
D.P = zeros(3,3,ne);  D.V = zeros(3,3,ne);  D.gm = zeros(ne,1);
D.mT = zeros(2,ne);   D.mF = zeros(2,ne);
for e = 1:ne
    t = datetime(char(ue(e)),'InputFormat','dd MMM yyyy HH:mm:ss.SSS','Locale','en_US');
    jd = datenum(t) + 1721058.5;
    D.gm(e) = deg2rad(mod(280.46061837 + 360.98564736629*(jd-2451545.0), 360));
    for k = 1:3
        i = find(S.Epoch==ue(e) & S.Satellite==sats(k), 1);
        if isempty(i), error('Missing state: %s at %s', sats(k), ue(e)); end
        D.P(k,:,e) = [S.X(i) S.Y(i) S.Z(i)];
        D.V(k,:,e) = [S.VX(i) S.VY(i) S.VZ(i)];
    end
    for k = 2:3
        j = find(M.Epoch==ue(e) & M.Satellite==sats(k), 1);
        if isempty(j), error('Missing measurement: %s at %s', sats(k), ue(e)); end
        D.mT(k-1,e) = M.Measured_TDOA_s(j);
        D.mF(k-1,e) = M.Measured_FDOA_Hz(j);
    end
end
fprintf('\n===== MODULE 09: CUMULATIVE SOLVER =====\n');
fprintf('Epochs: %d | measurements per epoch: 4 (2 TDOA + 2 FDOA)\n', ne);

%% 2. Search-grid centre from the satellite sub-points at epoch 1
g = D.gm(1);  Rz = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
lat0 = zeros(1,3); lon0 = zeros(1,3);
for k = 1:3
    pe = (Rz' * D.P(k,:,1)')';                 % inertial -> ECEF
    lat0(k) = asind(pe(3)/norm(pe));  lon0(k) = atan2d(pe(2),pe(1));
end
ctr = [mean(lat0) mean(lon0)];
fprintf('Grid centre (sub-point mean): lat %.2f, lon %.2f\n', ctr);
[gLat,gLon] = meshgrid(ctr(1)+(-cfg.gridHalf:cfg.gridStep:cfg.gridHalf), ...
                       ctr(2)+(-cfg.gridHalf:cfg.gridStep:cfg.gridHalf));
gLat = gLat(:);  gLon = gLon(:);
useLSQ = license('test','Optimization_Toolbox') && exist('lsqnonlin','file');
if ~useLSQ, fprintf('Optimization Toolbox not found: using fminsearch.\n'); end
opt = optimoptions_safe(useLSQ);

%% 3. Cumulative solve: epochs 1..n
rows = zeros(ne, 11);  xPrev = [];
for n = 1:ne
    ep = 1:n;
    cost = zeros(numel(gLat),1);
    for q = 1:numel(gLat)
        r = resid([gLat(q) gLon(q)], ep, D, cfg);
        cost(q) = sum(r.^2);
    end
    [~,ord] = sort(cost);
    starts = [gLat(ord(1:3)) gLon(ord(1:3))];
    if ~isempty(xPrev), starts = [starts; xPrev]; end

    best = Inf;  xBest = starts(1,:);
    for s = 1:size(starts,1)
        f = @(x) resid(x, ep, D, cfg);
        if useLSQ
            x = lsqnonlin(f, starts(s,:), [-90 -180], [90 180], opt);
        else
            x = fminsearch(@(x) sum(f(x).^2), starts(s,:), opt);
        end
        cst = sum(f(x).^2);
        if cst < best, best = cst; xBest = x; end
    end
    xPrev = xBest;

    % Covariance from numeric Jacobian (central differences)
    h = 1e-5;  J = zeros(4*n, 2);
    for p = 1:2
        dx = zeros(1,2); dx(p) = h;
        J(:,p) = (resid(xBest+dx,ep,D,cfg) - resid(xBest-dx,ep,D,cfg))/(2*h);
    end
    JtJ = J'*J;
    if rcond(JtJ) > 1e-14
        C = inv(JtJ);                                   % deg^2
        kmLat = 111.19;  kmLon = 111.19*cosd(xBest(1));
        sLat = sqrt(C(1,1))*kmLat;  sLon = sqrt(C(2,2))*kmLon;
    else
        sLat = NaN; sLon = NaN;
    end

    % Error vs truth (scoring only)
    eT = geo2ecef(cfg.truth(1),cfg.truth(2),cfg.truth(3));
    eE = geo2ecef(xBest(1),xBest(2),cfg.assumedAlt_km);
    errKm = norm(eE - eT);
    dof = 4*n - 2;
    rows(n,:) = [n 4*n xBest best dof best/max(dof,1) sLat sLon errKm cond(JtJ)];
    fprintf('n=%d | lat %.4f lon %.4f | chi2/dof %.2f | 1-sigma (lat,lon) %.2f, %.2f km | error %.2f km\n', ...
        n, xBest(1), xBest(2), best/max(dof,1), sLat, sLon, errKm);
end

%% 4. Save and summarize
names = {'NumEpochs','NumMeas','Lat_deg','Lon_deg','Chi2','DoF','Chi2_per_DoF', ...
         'SigmaLat_km','SigmaLon_km','Error_km','Cond_JtJ'};
R = array2table(rows(:,[1 2 3 4 5 6 7 8 9 10 11]), 'VariableNames', names);
writetable(R, 'Module09_LocalizationResults.csv');
fprintf('\nFinal estimate (n=%d): lat %.4f, lon %.4f | truth %.2f, %.2f | error %.2f km\n', ...
    ne, rows(end,3), rows(end,4), cfg.truth(1), cfg.truth(2), rows(end,10));
fprintf('Saved: Module09_LocalizationResults.csv\nMODULE 09 FINISHED.\n');

%% Local functions (keep at end of script)
function r = resid(x, ep, D, cfg)
    c = cfg.c;
    rE = geo2ecef(x(1), x(2), cfg.assumedAlt_km);
    r = zeros(4*numel(ep),1);
    for q = 1:numel(ep)
        e = ep(q);  g = D.gm(e);
        Rz = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
        re = (Rz*rE')';
        ve = cross([0 0 cfg.omegaE], re);
        rho = zeros(1,3);  rr = zeros(1,3);
        for k = 1:3
            dr = re - D.P(k,:,e);  dv = ve - D.V(k,:,e);
            rho(k) = norm(dr);  rr(k) = dot(dr,dv)/rho(k);
        end
        T = (rho(2:3) - rho(1))/c;
        F = cfg.fdoaSign*cfg.fc_Hz/c*(rr(2:3) - rr(1));
        r(4*q-3:4*q) = [(T' - D.mT(:,e))/cfg.sigT_s; (F' - D.mF(:,e))/cfg.sigF_Hz];
    end
end

function p = geo2ecef(latDeg, lonDeg, alt)
    a = 6378.137; fl = 1/298.257223563; e2 = fl*(2-fl);
    phi = deg2rad(latDeg); lam = deg2rad(lonDeg);
    N = a/sqrt(1 - e2*sin(phi)^2);
    p = [(N+alt)*cos(phi)*cos(lam), (N+alt)*cos(phi)*sin(lam), (N*(1-e2)+alt)*sin(phi)];
end

function opt = optimoptions_safe(useLSQ)
    if useLSQ
        opt = optimoptions('lsqnonlin','Display','off','TolFun',1e-14, ...
            'TolX',1e-12,'MaxIterations',400);
    else
        opt = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',2000);
    end
end