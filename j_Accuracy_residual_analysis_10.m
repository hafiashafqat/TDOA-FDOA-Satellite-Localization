%% ========================================================================
%% MODULE 10: ACCURACY, RESIDUAL, CONVERGENCE AND INITIAL-GUESS ANALYSIS
%% ========================================================================
% Needs: Module05_SatelliteStates.csv, Module07_SyntheticMeasurements.csv,
%        Module09_LocalizationResults.csv (optional, for the error-vs-n plot)
% Output: Module10_StartSensitivity.csv, Module10_ResidualSummary.csv
clc; format short g;

%% CONFIG (must match Modules 6, 7 and 9)
cfg.m5file = 'Module05_SatelliteStates.csv';
cfg.m7file = 'Module07_SyntheticMeasurements.csv';
cfg.m9file = 'Module09_LocalizationResults.csv';
cfg.fc_Hz = 2.4e9;  cfg.sigT_s = 1e-6;  cfg.sigF_Hz = 1.0;
cfg.fdoaSign = -1;  cfg.assumedAlt_km = 0.2;  cfg.omegaE = 7.2921150e-5;
cfg.truth = [31.00 73.50 0.2];            % scoring ONLY
cfg.successKm = 1.0;                      % a start "succeeds" if error < this
cfg.nRandomStarts = 20;
c = 299792.458;  cfg.c = c;  sats = ["AutoSat1","AutoSat2","AutoSat3"];
fprintf('\n===== MODULE 10: ACCURACY / RESIDUAL / CONVERGENCE =====\n');

%% 1. Load data (same as Module 9)
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
        D.P(k,:,e) = [S.X(i) S.Y(i) S.Z(i)];  D.V(k,:,e) = [S.VX(i) S.VY(i) S.VZ(i)];
    end
    for k = 2:3
        j = find(M.Epoch==ue(e) & M.Satellite==sats(k), 1);
        D.mT(k-1,e) = M.Measured_TDOA_s(j);  D.mF(k-1,e) = M.Measured_FDOA_Hz(j);
    end
end
useLSQ = license('test','Optimization_Toolbox') && exist('lsqnonlin','file');
opt = optimoptions_safe(useLSQ);
errFun = @(x) norm(geo2ecef(x(1),x(2),cfg.assumedAlt_km) - ...
                   geo2ecef(cfg.truth(1),cfg.truth(2),cfg.truth(3)));

% Grid centre from epoch-1 sub-points (same rule as Module 9)
g = D.gm(1);  Rz = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
la = zeros(1,3); lo = zeros(1,3);
for k = 1:3
    pe = (Rz'*D.P(k,:,1)')';  la(k) = asind(pe(3)/norm(pe));  lo(k) = atan2d(pe(2),pe(1));
end
ctr = [mean(la) mean(lo)];

%% 2. Reference solution on the full arc (grid search + refinement)
ep = 1:ne;
[gLat,gLon] = meshgrid(ctr(1)+(-12:1:12), ctr(2)+(-12:1:12));
gLat = gLat(:); gLon = gLon(:);
cg = arrayfun(@(a,b) sum(resid([a b],ep,D,cfg).^2), gLat, gLon);
[~,ord] = sort(cg);  bestC = Inf;  xHat = [ctr];
for s = 1:3
    [x,cst] = solveFrom([gLat(ord(s)) gLon(ord(s))], ep, D, cfg, useLSQ, opt);
    if cst < bestC, bestC = cst; xHat = x; end
end

%% 3. Accuracy at the full arc
h = 1e-5; J = zeros(4*ne,2);
for p = 1:2
    dx = zeros(1,2); dx(p) = h;
    J(:,p) = (resid(xHat+dx,ep,D,cfg) - resid(xHat-dx,ep,D,cfg))/(2*h);
end
C = inv(J'*J);                                        % deg^2
kmLat = 111.19;  kmLon = 111.19*cosd(xHat(1));
eN = (xHat(1)-cfg.truth(1))*kmLat;  eE = (xHat(2)-cfg.truth(2))*kmLon;
eVec = [xHat(1)-cfg.truth(1); xHat(2)-cfg.truth(2)];
d2 = eVec'*(J'*J)*eVec;                               % expect ~2 (chi2, 2 dof)
fprintf('\n--- ACCURACY (n=%d) ---\n', ne);
fprintf('Estimate: lat %.5f, lon %.5f | truth %.2f, %.2f\n', xHat, cfg.truth(1:2));
fprintf('Error: north %.4f km, east %.4f km, total %.4f km\n', eN, eE, errFun(xHat));
fprintf('Formal 1-sigma: north %.4f km, east %.4f km\n', sqrt(C(1,1))*kmLat, sqrt(C(2,2))*kmLon);
fprintf('Mahalanobis^2 of the error: %.2f (expect ~2; single realization, so only indicative)\n', d2);
fprintf('Jacobian condition number: %.1f\n', cond(J'*J));

%% 4. Residual analysis at the estimate
r  = resid(xHat, ep, D, cfg);   R4 = reshape(r,4,ne)';   % [T21 T31 F21 F31] / sigma
rT = resid([cfg.truth(1) cfg.truth(2)], ep, D, cfg);
dof = 4*ne - 2;  chi2 = sum(r.^2);  chi2T = sum(rT.^2);
if exist('chi2inv','file')
    lo95 = chi2inv(0.025,dof);  hi95 = chi2inv(0.975,dof);
else
    lo95 = dof-1.96*sqrt(2*dof);  hi95 = dof+1.96*sqrt(2*dof);
end
lbl = {'TDOA 2-1','TDOA 3-1','FDOA 2-1','FDOA 3-1'};
unit = [c*1000*cfg.sigT_s, c*1000*cfg.sigT_s, cfg.sigF_Hz, cfg.sigF_Hz];  % m, m, Hz, Hz
fprintf('\n--- RESIDUALS ---\n');
fprintf('chi2 = %.2f, dof = %d, chi2/dof = %.2f (95%% band for chi2: %.1f to %.1f)\n', chi2, dof, chi2/dof, lo95, hi95);
fprintf('chi2 at truth = %.2f; at estimate = %.2f (estimate must not be larger)\n', chi2T, chi2);
fprintf('%-9s %9s %9s %12s %9s\n','Type','mean','std','max|res|','lag1');
sumRows = cell(4,5);
for q = 1:4
    x = R4(:,q);
    if ne > 2, ac = corr(x(1:end-1),x(2:end)); else, ac = NaN; end
    mu1 = mean(x); sd = std(x);
    if q <= 2, ustr = 'm'; else, ustr = 'Hz'; end
    fprintf('%-9s %9.3f %9.3f %8.3f (%.2f %s) %7.2f\n', lbl{q}, mu1, sd, max(abs(x)), max(abs(x))*unit(q), ustr, ac);
    sumRows(q,:) = {lbl{q}, mu1, sd, max(abs(x)), ac};
end
fprintf('(Normalized residuals should have mean ~0, std ~1, and small lag-1 correlation.)\n');
writetable(cell2table(sumRows,'VariableNames',{'Type','Mean','Std','MaxAbs','Lag1'}), ...
    'Module10_ResidualSummary.csv');

%% 5. Convergence and sensitivity to the initial guess
[sLat,sLon] = meshgrid(ctr(1)+(-15:5:15), ctr(2)+(-15:5:15));
rng(1);
rndLat = -60 + 120*rand(cfg.nRandomStarts,1);
rndLon = -180 + 360*rand(cfg.nRandomStarts,1);
starts = [sLat(:) sLon(:); rndLat rndLon];
starts(:,1) = min(max(starts(:,1),-89),89);
nS = size(starts,1);
caseN = [1 ne];  out = [];
fprintf('\n--- INITIAL-GUESS SENSITIVITY (%d starts per case, local solver only) ---\n', nS);
for cN = 1:numel(caseN)
    n = caseN(cN);  epn = 1:n;
    X = zeros(nS,2); Cst = zeros(nS,1); It = zeros(nS,1); Ef = zeros(nS,1); Er = zeros(nS,1);
    for s = 1:nS
        try
            [X(s,:),Cst(s),It(s),Ef(s)] = solveFrom(starts(s,:), epn, D, cfg, useLSQ, opt);
            Er(s) = errFun(X(s,:));
        catch
            X(s,:) = NaN; Cst(s) = Inf; It(s) = NaN; Ef(s) = NaN; Er(s) = Inf;
        end
    end
    ok = Er < cfg.successKm;
    fprintf('n=%d: success %d/%d (%.0f%%) | iterations median %.0f, max %.0f | cost of successes %.1f..%.1f\n', ...
        n, sum(ok), nS, 100*mean(ok), median(It,'omitnan'), max(It), ...
        min(Cst(ok)), max(Cst(ok)));
    if any(~ok)
        bad = find(~ok);
        fprintf('   %d failed starts. Local minima found (lat, lon, cost):\n', numel(bad));
        [~,u] = unique(round(X(bad,:),1),'rows','stable');
        for i = u(:)'
            fprintf('     start (%.1f, %.1f) -> (%.3f, %.3f), cost %.1f, error %.0f km\n', ...
                starts(bad(i),:), X(bad(i),:), Cst(bad(i)), Er(bad(i)));
        end
    end
    out = [out; table(repmat(n,nS,1), starts(:,1), starts(:,2), X(:,1), X(:,2), Cst, It, Ef, Er, ok, ...
        'VariableNames',{'NumEpochs','StartLat','StartLon','EndLat','EndLon','Cost','Iterations','ExitFlag','Error_km','Success'})]; %#ok<AGROW>
end
writetable(out,'Module10_StartSensitivity.csv');

%% 6. Plots
figure('Name','Module 10','Position',[80 80 1100 800]);
subplot(2,2,1);
if isfile(cfg.m9file)
    T9 = readtable(cfg.m9file,'VariableNamingRule','preserve');
    semilogy(T9.NumEpochs, max(T9.Error_km,1e-4),'o-','LineWidth',1.5); hold on;
    semilogy(T9.NumEpochs, hypot(T9.SigmaLat_km,T9.SigmaLon_km),'s--','LineWidth',1.5);
    legend('Actual error','Formal 1-sigma (radial)','Location','best');
end
xlabel('Epochs used'); ylabel('km'); title('Accuracy vs arc length'); grid on;
subplot(2,2,2);
plot(1:ne,R4,'o-'); yline(0,'k-'); yline([-3 3],'r:');
xlabel('Epoch'); ylabel('Normalized residual'); title('Residuals at the estimate'); legend(lbl,'Location','best'); grid on;
subplot(2,2,3);
for cN = 1:numel(caseN)
    sel = out.NumEpochs==caseN(cN);
    histogram(out.Iterations(sel),'DisplayStyle','stairs','LineWidth',1.5); hold on;
end
xlabel('Iterations to converge'); ylabel('Count'); title('Convergence'); legend('n=1',sprintf('n=%d',ne)); grid on;
subplot(2,2,4);
sel = out.NumEpochs==ne;
scatter(out.StartLon(sel),out.StartLat(sel),40,double(out.Success(sel)),'filled'); hold on;
plot(cfg.truth(2),cfg.truth(1),'kp','MarkerSize',14,'MarkerFaceColor','y');
colormap(subplot(2,2,4),[0.85 0.2 0.2; 0.2 0.7 0.3]);
xlabel('Start longitude'); ylabel('Start latitude'); title(sprintf('Start success map (n=%d); green = converged',ne)); grid on;

%% 7. Verdict
pRes = chi2 >= lo95 && chi2 <= hi95 && chi2 <= chi2T + 1e-6;
fprintf('\n===== VERDICT =====\n');
fprintf('Residual statistics consistent with the noise model: %s\n', ternary(pRes,'YES','NO - check model/sigma'));
fprintf('Total error %.3f km vs formal radial 1-sigma %.3f km\n', errFun(xHat), hypot(sqrt(C(1,1))*kmLat, sqrt(C(2,2))*kmLon));
fprintf('MODULE 10 FINISHED.\n');

%% Local functions (keep at end of script)
function [x,cst,it,ef] = solveFrom(x0, ep, D, cfg, useLSQ, opt)
    f = @(x) resid(x, ep, D, cfg);
    if useLSQ
        [x,cst,~,ef,o] = lsqnonlin(f, x0, [-90 -180], [90 180], opt);  it = o.iterations;
    else
        [x,cst,ef,o] = fminsearch(@(x) sum(f(x).^2), x0, opt);  it = o.iterations;
    end
end
function r = resid(x, ep, D, cfg)
    c = cfg.c;  rE = geo2ecef(x(1), x(2), cfg.assumedAlt_km);
    r = zeros(4*numel(ep),1);
    for q = 1:numel(ep)
        e = ep(q);  g = D.gm(e);
        Rz = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
        re = (Rz*rE')';  ve = cross([0 0 cfg.omegaE], re);
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
        opt = optimoptions('lsqnonlin','Display','off','TolFun',1e-14,'TolX',1e-12,'MaxIterations',400);
    else
        opt = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',2000);
    end
end
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end