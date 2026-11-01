%% ========================================================================
%% MODULE 11: MONTE CARLO AND ROBUSTNESS TESTING
%% ========================================================================
% Needs: Module05_SatelliteStates.csv
% Output: Module11_MonteCarlo.csv, Module11_Stress.csv
% Runtime: roughly 5-10 min at 100 trials. Set cfg.nTrials = 20 for a test run.
clc;
format short g;

%% CONFIG (must match Modules 6 and 7)
cfg.m5file   = 'Module05_SatelliteStates.csv';
cfg.fc_Hz    = 2.4e9;  cfg.fdoaSign = -1;  cfg.omegaE = 7.2921150e-5;
cfg.c        = 299792.458;
cfg.truthAlt = 0.2;                  % km, true emitter altitude (Module 6)
cfg.nTrials  = 100;                  % trials per case
cfg.seed     = 1000;                 % trial t uses rng(seed + t)
cfg.failKm   = 50;                   % error above this = gross failure
cfg.gridHalf = 25;  cfg.gridStep = 5;  cfg.nStarts = 3;   % search grid (deg)
sats = ["AutoSat1","AutoSat2","AutoSat3"];

%% 1. Load satellite states
if ~isfile(cfg.m5file), error('%s not found in current folder.', cfg.m5file); end
S = readtable(cfg.m5file, setvartype(detectImportOptions(cfg.m5file, ...
    'VariableNamingRule','preserve'),'Epoch','string'));
S.Epoch = strtrim(string(S.Epoch));  S.Satellite = strtrim(string(S.Satellite));
ue = unique(S.Epoch,'stable');  ne = numel(ue);
D.P = zeros(3,3,ne);  D.V = zeros(3,3,ne);  D.gm = zeros(ne,1);
for e = 1:ne
    tt = datetime(char(ue(e)),'InputFormat','dd MMM yyyy HH:mm:ss.SSS','Locale','en_US');
    jd = datenum(tt) + 1721058.5;
    D.gm(e) = deg2rad(mod(280.46061837 + 360.98564736629*(jd-2451545.0), 360));
    for k = 1:3
        ii = find(S.Epoch==ue(e) & S.Satellite==sats(k), 1);
        if isempty(ii), error('Missing state: %s at %s', sats(k), ue(e)); end
        D.P(k,:,e) = [S.X(ii) S.Y(ii) S.Z(ii)];
        D.V(k,:,e) = [S.VX(ii) S.VY(ii) S.VZ(ii)];
    end
end
cfg.useLSQ = license('test','Optimization_Toolbox') && exist('lsqnonlin','file');
cfg.opt = optimoptions_safe(cfg.useLSQ);
if ~cfg.useLSQ, fprintf('Optimization Toolbox not found: using fminsearch (slower).\n'); end

% Search-grid centre from epoch-1 sub-points (same rule as Module 9)
g1 = D.gm(1);  Rz1 = [cos(g1) -sin(g1) 0; sin(g1) cos(g1) 0; 0 0 1];
subLat = zeros(1,3);  subLon = zeros(1,3);
for k = 1:3
    pe = (Rz1'*D.P(k,:,1)')';
    subLat(k) = asind(pe(3)/norm(pe));  subLon(k) = atan2d(pe(2),pe(1));
end
ctr = [mean(subLat) mean(subLon)];

%% 2. Define geometry cases and noise levels
epShort  = 1:min(3,ne);
epSparse = unique([1 round(ne/2) ne]);
geoms = { ...
    'Baseline (ref 1, all epochs)', [1 2 3], 1:ne,     [31.0 73.5]; ...
    'Reference = AutoSat2',         [2 1 3], 1:ne,     [31.0 73.5]; ...
    'Reference = AutoSat3',         [3 1 2], 1:ne,     [31.0 73.5]; ...
    'Short arc (first 3 epochs)',   [1 2 3], epShort,  [31.0 73.5]; ...
    'Sparse arc (start/mid/end)',   [1 2 3], epSparse, [31.0 73.5]; ...
    'Two satellites (1,2)',         [1 2],   1:ne,     [31.0 73.5]; ...
    'Emitter 2 deg north',          [1 2 3], 1:ne,     [33.0 73.5]; ...
    'Emitter 2 deg south',          [1 2 3], 1:ne,     [29.0 73.5]; ...
    'Emitter 2.5 deg east',         [1 2 3], 1:ne,     [31.0 76.0]};
nG = size(geoms,1);
noise = struct('name',{'0.1us/0.1Hz','1us/1Hz','10us/10Hz','100us/100Hz'}, ...
               'sigT',{1e-7,1e-6,1e-5,1e-4}, 'sigF',{0.1,1,10,100});
nN = numel(noise);
noStress = struct('altOff',0,'fcErr',0);

fprintf('\n===== MODULE 11: MONTE CARLO (%d trials per case) =====\n', cfg.nTrials);
fprintf('Epochs available: %d | search grid centre: lat %.2f, lon %.2f\n', ne, ctr);

%% 3. Noise x geometry sweep
rows = {};  rmsM = nan(nG,nN);  ratM = nan(nG,nN);
tStart = tic;  firstDone = false;
for gi = 1:nG
    gs = struct('name',geoms{gi,1}, 'si',geoms{gi,2}, 'ep',geoms{gi,3}, ...
                'truth',[geoms{gi,4} cfg.truthAlt]);
    el = minElev(gs, D);
    if el < 0
        fprintf('\n[SKIP] %s: a satellite is below the horizon (min elevation %.1f deg)\n', gs.name, el);
        continue;
    end
    fprintf('\n%s | min elevation %.1f deg | %d sats, %d epochs\n', ...
        gs.name, el, numel(gs.si), numel(gs.ep));
    for ni = 1:nN
        res = runCell(gs, noise(ni), noStress, D, cfg, ctr, cfg.nTrials);
        sm  = summarize(res, cfg);
        rmsM(gi,ni) = sm.rms;  ratM(gi,ni) = sm.ratio;
        fprintf('  %-12s fail %5.1f%% | median %9.4f km | RMS %9.4f | P95 %9.4f | emp/formal %5.2f | cov95 %3.0f%% | chi2/dof %.2f\n', ...
            noise(ni).name, sm.failPct, sm.median, sm.rms, sm.p95, sm.ratio, sm.cov95, sm.chi2);
        rows(end+1,:) = {gs.name, noise(ni).name, el, cfg.nTrials, sm.failPct, sm.median, ...
            sm.rms, sm.p95, sm.formal, sm.ratio, sm.cov95, sm.chi2}; %#ok<AGROW>
        if ~firstDone
            firstDone = true;
            fprintf('  (first case took %.0f s; estimated total about %.0f min)\n', ...
                toc(tStart), toc(tStart)*(nG*nN + 8)/60);
        end
    end
end
if isempty(rows), error('All geometry cases were skipped.'); end
T11 = cell2table(rows,'VariableNames',{'Geometry','Noise','MinElev_deg','Trials','FailPct', ...
    'Median_km','RMS_km','P95_km','FormalRMS_km','EmpOverFormal','Cov95Pct','Chi2PerDoF'});
writetable(T11,'Module11_MonteCarlo.csv');

%% 4. Stress tests on the baseline case (baseline noise, 1 us / 1 Hz)
gB = struct('name',geoms{1,1}, 'si',geoms{1,2}, 'ep',geoms{1,3}, ...
            'truth',[geoms{1,4} cfg.truthAlt]);
altOff = [0 0.5 1 3];  fcErr = [0 1e-6 1e-5 1e-4];
sRows = {};  sRmsAlt = zeros(1,numel(altOff));  sRmsFc = zeros(1,numel(fcErr));
fprintf('\n--- STRESS: assumed-altitude error (truth %.1f km) ---\n', cfg.truthAlt);
for i = 1:numel(altOff)
    res = runCell(gB, noise(2), struct('altOff',altOff(i),'fcErr',0), D, cfg, ctr, cfg.nTrials);
    sm = summarize(res, cfg);  sRmsAlt(i) = sm.rms;
    fprintf('  alt error %+4.1f km | fail %5.1f%% | median %8.4f km | RMS %8.4f | cov95 %3.0f%% | chi2/dof %.2f\n', ...
        altOff(i), sm.failPct, sm.median, sm.rms, sm.cov95, sm.chi2);
    sRows(end+1,:) = {'AltitudeError_km', altOff(i), sm.failPct, sm.median, sm.rms, sm.cov95, sm.chi2}; %#ok<AGROW>
end
fprintf('\n--- STRESS: carrier-frequency relative error (generation vs solver) ---\n');
for i = 1:numel(fcErr)
    res = runCell(gB, noise(2), struct('altOff',0,'fcErr',fcErr(i)), D, cfg, ctr, cfg.nTrials);
    sm = summarize(res, cfg);  sRmsFc(i) = sm.rms;
    fprintf('  fc error %.0e (%.1f kHz) | fail %5.1f%% | median %8.4f km | RMS %8.4f | cov95 %3.0f%% | chi2/dof %.2f\n', ...
        fcErr(i), fcErr(i)*cfg.fc_Hz/1e3, sm.failPct, sm.median, sm.rms, sm.cov95, sm.chi2);
    sRows(end+1,:) = {'CarrierRelError', fcErr(i), sm.failPct, sm.median, sm.rms, sm.cov95, sm.chi2}; %#ok<AGROW>
end
writetable(cell2table(sRows,'VariableNames', ...
    {'Stress','Value','FailPct','Median_km','RMS_km','Cov95Pct','Chi2PerDoF'}), 'Module11_Stress.csv');

%% 5. Plots
figure('Name','Module 11','Position',[80 80 1150 800]);
xs = [noise.sigT]*1e6;
subplot(2,2,1);
loglog(xs, rmsM','o-','LineWidth',1.3); grid on;
xlabel('TDOA sigma (us); FDOA sigma scales with it'); ylabel('RMS horizontal error (km)');
title('Error vs noise level'); legend(geoms(:,1),'Location','southeast','FontSize',7);
subplot(2,2,2);
semilogx(xs, ratM','o-','LineWidth',1.3); yline(1,'k--'); grid on;
xlabel('TDOA sigma (us)'); ylabel('Empirical RMS / formal 1-sigma');
title('Covariance calibration (1 = ideal)');
subplot(2,2,3);
bar(categorical(string(altOff)), sRmsAlt); set(gca,'YScale','log'); grid on;
xlabel('Assumed altitude error (km)'); ylabel('RMS error (km)'); title('Altitude error sensitivity');
subplot(2,2,4);
bar(categorical(string(fcErr)), sRmsFc); set(gca,'YScale','log'); grid on;
xlabel('Carrier relative error'); ylabel('RMS error (km)'); title('Carrier-frequency error sensitivity');

%% 6. Verdict
fprintf('\n===== VERDICT =====\n');
sel = strcmp(T11.Geometry, geoms{1,1}) & strcmp(T11.Noise, '1us/1Hz');
if any(sel)
    b = T11(find(sel,1),:);
    fprintf('Baseline at Module 7 noise: fail %.1f%%, RMS %.4f km, empirical/formal %.2f, cov95 %.0f%%, chi2/dof %.2f\n', ...
        b.FailPct, b.RMS_km, b.EmpOverFormal, b.Cov95Pct, b.Chi2PerDoF);
    okCal = b.EmpOverFormal > 0.8 && b.EmpOverFormal < 1.25 && ...
            b.Chi2PerDoF > 0.8 && b.Chi2PerDoF < 1.2;
    if okCal, msg = 'YES'; else, msg = 'NO - inspect'; end
    fprintf('Formal covariance calibrated against Monte Carlo scatter: %s\n', msg);
end
fprintf('Total runtime: %.1f min\nMODULE 11 FINISHED.\n', toc(tStart)/60);

%% Local functions (must stay at the end of the script)
function res = runCell(gs, nz, st, D, cfg, ctr, nT)
    alt = cfg.truthAlt;  ep = gs.ep;  si = gs.si;  m = numel(si)-1;  nEp = numel(ep);
    aAlt = alt + st.altOff;  fcTrue = cfg.fc_Hz*(1 + st.fcErr);
    [T0,F0] = predict(gs.truth(1), gs.truth(2), alt, ep, si, D, cfg, fcTrue);
    [gLat,gLon] = meshgrid(ctr(1)+(-cfg.gridHalf:cfg.gridStep:cfg.gridHalf), ...
                           ctr(2)+(-cfg.gridHalf:cfg.gridStep:cfg.gridHalf));
    gLat = min(max(gLat(:),-89),89);  gLon = gLon(:);
    dof = 2*m*nEp - 2;
    tr  = geo2ecef(gs.truth(1), gs.truth(2), alt);
    res.err = nan(nT,1);  res.frm = nan(nT,1);  res.d2 = nan(nT,1);  res.c2 = nan(nT,1);
    cs.ep = ep;  cs.si = si;  cs.alt = aAlt;  cs.sigT = nz.sigT;  cs.sigF = nz.sigF;
    for t = 1:nT
        rng(cfg.seed + t);
        Zt = randn(2,nEp);  Zf = randn(2,nEp);            % same draws in every case
        cs.mT = T0 + nz.sigT*Zt(1:m,:);
        cs.mF = F0 + nz.sigF*Zf(1:m,:);
        cost = arrayfun(@(a,b) sum(residual([a b],cs,D,cfg).^2), gLat, gLon);
        [~,ord] = sort(cost);
        best = Inf;  xb = [NaN NaN];
        for s = 1:cfg.nStarts
            try
                [x,cst] = solveFrom([gLat(ord(s)) gLon(ord(s))], cs, D, cfg);
            catch
                continue;
            end
            if cst < best, best = cst; xb = x; end
        end
        if ~isfinite(best), continue; end
        res.err(t) = norm(geo2ecef(xb(1), xb(2), alt) - tr);   % horizontal error at true altitude
        res.c2(t)  = best/max(dof,1);
        h = 1e-5;  J = zeros(numel(residual(xb,cs,D,cfg)),2);
        for p = 1:2
            dx = zeros(1,2);  dx(p) = h;
            J(:,p) = (residual(xb+dx,cs,D,cfg) - residual(xb-dx,cs,D,cfg))/(2*h);
        end
        JtJ = J'*J;
        if rcond(JtJ) > 1e-14
            C = inv(JtJ);                                      % deg^2
            res.frm(t) = hypot(sqrt(C(1,1))*111.19, sqrt(C(2,2))*111.19*cosd(xb(1)));
            ev = [xb(1)-gs.truth(1); xb(2)-gs.truth(2)];
            res.d2(t) = ev'*JtJ*ev;                            % chi2, 2 dof
        end
    end
end

function sm = summarize(res, cfg)
    fail = ~(res.err < cfg.failKm);   ok = ~fail;   e = res.err(ok);
    sm.failPct = 100*mean(fail);
    if isempty(e)
        sm.median = NaN; sm.rms = NaN; sm.p95 = NaN; sm.formal = NaN;
        sm.ratio = NaN; sm.cov95 = NaN; sm.chi2 = NaN;
        return;
    end
    s = sort(e);
    sm.median = median(e);  sm.rms = sqrt(mean(e.^2));
    sm.p95 = s(max(1,ceil(0.95*numel(s))));
    f = res.frm(ok);  f = f(~isnan(f));
    sm.formal = sqrt(mean(f.^2));  sm.ratio = sm.rms/sm.formal;
    d = res.d2(ok);  d = d(~isnan(d));
    sm.cov95 = 100*mean(d < 5.991);                            % chi2inv(0.95,2)
    sm.chi2  = mean(res.c2(ok),'omitnan');
end

function el = minElev(gs, D)
    rE = geo2ecef(gs.truth(1), gs.truth(2), gs.truth(3));  el = Inf;
    for e = gs.ep
        gm = D.gm(e);  Rz = [cos(gm) -sin(gm) 0; sin(gm) cos(gm) 0; 0 0 1];
        re = (Rz*rE')';
        for k = gs.si
            dr = D.P(k,:,e) - re;
            el = min(el, asind(dot(dr, re/norm(re))/norm(dr)));
        end
    end
end

function [T,F] = predict(lat, lon, alt, ep, si, D, cfg, fc)
    c = cfg.c;  rE = geo2ecef(lat, lon, alt);  m = numel(si)-1;
    T = zeros(m,numel(ep));  F = T;
    for q = 1:numel(ep)
        e = ep(q);  gm = D.gm(e);
        Rz = [cos(gm) -sin(gm) 0; sin(gm) cos(gm) 0; 0 0 1];
        re = (Rz*rE')';  ve = cross([0 0 cfg.omegaE], re);
        rho = zeros(1,m+1);  rr = zeros(1,m+1);
        for k = 1:m+1
            dr = re - D.P(si(k),:,e);  dv = ve - D.V(si(k),:,e);
            rho(k) = norm(dr);  rr(k) = dot(dr,dv)/rho(k);
        end
        T(:,q) = ((rho(2:end) - rho(1))/c)';
        F(:,q) = (cfg.fdoaSign*fc/c*(rr(2:end) - rr(1)))';
    end
end

function r = residual(x, cs, D, cfg)
    [T,F] = predict(x(1), x(2), cs.alt, cs.ep, cs.si, D, cfg, cfg.fc_Hz);
    r = [(T(:) - cs.mT(:))/cs.sigT; (F(:) - cs.mF(:))/cs.sigF];
end

function [x,cst] = solveFrom(x0, cs, D, cfg)
    f = @(x) residual(x, cs, D, cfg);
    if cfg.useLSQ
        [x,cst] = lsqnonlin(f, x0, [-90 -180], [90 180], cfg.opt);
    else
        [x,cst] = fminsearch(@(x) sum(f(x).^2), x0, cfg.opt);
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
        opt = optimoptions('lsqnonlin','Display','off','TolFun',1e-14,'TolX',1e-12,'MaxIterations',400);
    else
        opt = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',2000);
    end
end