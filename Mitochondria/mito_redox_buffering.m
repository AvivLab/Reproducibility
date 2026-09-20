function mito_redox_buffering()
% Toy dynamical model INSPIRED BY the mitochondrial energy-redox chaos
% result reported in
%   Kembro, Cortassa, Lloyd, Sollott & Aon, 2018, "Mitochondrial chaotic
%   dynamics: Redox-energetic behavior at the edge of stability,"
%   Scientific Reports 8:15422.
% That paper's own model (the Cortassa/Aon "ME-R" mitochondrial
% energy-redox-Ca2+ model) has dozens of state variables: respiratory
% chain complexes I-IV, the TCA cycle, ANT, F1F0-ATPase, Ca2+ handling,
% GSH/GSSG and thioredoxin pools, and separate mitochondrial (SOD2) and
% cytoplasmic/intermembrane-space (SOD1) superoxide dismutase fluxes.
% Reproducing that model exactly is well beyond a toy illustration, and
% nothing here should be read as a reconstruction of it. What IS carried
% over faithfully is the qualitative finding the paper reports: reactive
% oxygen species (ROS) handling and mitochondrial bioenergetics are
% mutually coupled through cross-talk between the redox state and the
% respiratory/TCA machinery, antioxidant (SOD) buffering capacity
% ordinarily keeps that coupled system out of chaos, and the system sits
% "at the edge of stability" -- small reductions in buffering capacity
% can flip it into genuine deterministic chaos (the paper reports
% positive Lyapunov exponents and a strange attractor), not through a
% single clean on/off switch but through a sensitive, non-monotonic
% dependence on how much buffering capacity remains.
%
% Network (same two-ring cross-talk architecture used throughout this
% project's other toy models, here relabeled for the redox/bioenergetic
% story rather than re-derived from scratch):
%   ROS-handling ring:     Sox -| Prx -| Rdx -| Sox
%     (superoxide -| peroxide/peroxiredoxin flux -| overall thiol/
%      glutathione redox hub)
%   Bioenergetic ring:     Psi -| Atp -| Suc -| Psi
%     (membrane potential -| ATP/ADP ratio -| TCA/succinate flux hub)
%   Cross-talk (real coupling; oxidative modification of TCA enzymes and
%     Complex II/reverse-electron-transport ROS production genuinely
%     link the two, sub-threshold on its own): Rdx <-> Suc
% Two antioxidant buffering fluxes, Sod1 (cytoplasmic/intermembrane
% -space-like) and Sod2 (mitochondrial-matrix-like), each activated by
% the redox hub and each adding a clearance term onto BOTH cross-talk
% hubs -- a fast repair/scavenging action that damps runaway excursions,
% not an extra destabilizing signal. Their COMBINED strength is the
% single "buffering capacity" parameter (k) that this script varies.
%
% A parameter scan (see buf_sensitivity.m in this same folder) confirmed
% the network is intrinsically chaotic when k = 0 (bare cross-talk,
% c = 3.0, same chaos threshold established for the other two-ring toy
% models in this project) and mapped out how chaos depends on k across
% k = 0 to 0.15: chaos persists (with one narrow k = 0.02 stability
% island) until k reaches roughly 0.06-0.07, above which the network is
% robustly stable for the rest of the range tested. This non-monotonic,
% threshold-like sensitivity -- rather than a smooth, gradual quenching
% of chaos as buffering increases -- is itself a reasonable toy analogue
% of the "edge of stability" behavior Kembro et al. describe, where
% "slight parametric changes... lead to drastic qualitative changes in
% dynamics."
%
% Two conditions compared below:
%   normal buffering   (k = 0.10, safely inside the stable region)
%   reduced buffering  (k = 0.01, deep in the chaotic region)
% Buffering is reduced, not eliminated, on purpose: even a small residual
% loss of capacity is enough to reveal the chaos, echoing the real
% paper's point that this is a sensitive edge, not a simple present/
% absent switch.

close all;

n  = 4;  K = 1;
gA = 1.0; gB = 2.0;
aA = 8.0; aB = 8.0;
c  = 3.0;                     % Rdx<->Suc cross-talk, above the chaos threshold
gBuf = 1.5; aBuf = 8.0;       % Sod1, Sod2: identical turnover and gain

Hact = @(x) (max(x,0).^n) ./ (K^n + max(x,0).^n);
Hrep = @(x) 1 - Hact(x);

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.5);
Tspan = [0 300];
tt = linspace(0, 300, 3000);

x0 = [0.42 0.72 0.10 0.68 0.24 0.19 0.30 0.55]';
x1 = x0; x1(1) = x1(1) + 1e-6;

k_normal  = 0.10;
k_reduced = 0.01;

rhs_normal  = @(t,x) buffered_rhs(t,x,gA,gB,aA,aB,c,gBuf,aBuf,k_normal,Hact,Hrep);
rhs_reduced = @(t,x) buffered_rhs(t,x,gA,gB,aA,aB,c,gBuf,aBuf,k_reduced,Hact,Hrep);

[t0n,y0n] = ode45(rhs_normal, Tspan, x0, opts);
[t1n,y1n] = ode45(rhs_normal, Tspan, x1, opts);
Y0n = interp1(t0n,y0n,tt); Y1n = interp1(t1n,y1n,tt);
sep_n = sqrt(sum((Y0n-Y1n).^2,2));

[t0r,y0r] = ode45(rhs_reduced, Tspan, x0, opts);
[t1r,y1r] = ode45(rhs_reduced, Tspan, x1, opts);
Y0r = interp1(t0r,y0r,tt); Y1r = interp1(t1r,y1r,tt);
sep_r = sqrt(sum((Y0r-Y1r).^2,2));

fprintf('NORMAL buffering (k=%.2f):  sep(end)=%.6f  max=%.6f\n', k_normal, sep_n(end), max(sep_n));
fprintf('REDUCED buffering (k=%.2f): sep(end)=%.6f  max=%.6f\n', k_reduced, sep_r(end), max(sep_r));

% ---------------------------------------------------------------
% Figure 1: buffering-strength sensitivity scan (from buf_sensitivity.m)
% ---------------------------------------------------------------
if exist('sensitivity_data.mat', 'file')
  load('sensitivity_data.mat');
  fig1 = figure('Position', [80 80 700 400]);
  semilogy(kvals, max(seps,1e-7), 'o-', 'Color', [0.55 0.05 0.05], 'LineWidth', 1.3, 'MarkerFaceColor', [0.55 0.05 0.05]);
  hold on;
  xlim([-0.01 0.16]);
  plot([k_normal k_normal], ylim, 'b--', 'LineWidth', 1.2);
  plot([k_reduced k_reduced], ylim, '--', 'Color', [0.85 0.45 0], 'LineWidth', 1.2);
  text(k_normal+0.003, 3, 'normal buffering', 'Color', 'b');
  text(k_reduced+0.003, 8e-7*100, 'reduced buffering', 'Color', [0.85 0.45 0]);
  xlabel('combined Sod1+Sod2 buffering strength, k');
  ylabel('trajectory separation at t=300');
  title('Chaos is a sensitive, non-monotonic function of buffering capacity');
  box on;
  print(fig1, 'mito_buffering_sensitivity.png', '-dpng', '-r150');
  fprintf('Figure saved to mito_buffering_sensitivity.png\n');
end

% ---------------------------------------------------------------
% Figure 2: before/after divergence test and hub phase portrait
% ---------------------------------------------------------------
fig2 = figure('Position', [100 100 950 400]);

subplot(1,2,1);
semilogy(tt, sep_n, 'b-', 'LineWidth', 1.4); hold on;
semilogy(tt, sep_r, '-', 'Color', [0.85 0.45 0], 'LineWidth', 1.4);
xlabel('time'); ylabel('separation between nearby trajectories');
legend('normal buffering (k=0.10)', 'reduced buffering (k=0.01)', 'Location', 'northwest');
title('Divergence test, before vs after buffering reduction');
box on;

tail = tt > 100;
subplot(1,2,2);
plot(Y0r(tail,3), Y0r(tail,6), '-', 'Color', [0.85 0.45 0], 'LineWidth', 0.5); hold on;
plot(Y0n(tail,3), Y0n(tail,6), 'b-', 'LineWidth', 1.4);
xlabel('Rdx (redox hub)'); ylabel('Suc (TCA/succinate hub)');
legend('reduced buffering (chaotic)', 'normal buffering (stable)', 'Location', 'best');
title('Hub phase portrait, t>100');
axis square; box on;

print(fig2, 'mito_buffering_divergence.png', '-dpng', '-r150');
fprintf('Figure saved to mito_buffering_divergence.png\n');

% ---------------------------------------------------------------
% Figure 3: Rdx(t) time series, before vs after, stacked
% ---------------------------------------------------------------
fig3 = figure('Position', [80 80 900 500]);
subplot(2,1,1);
plot(tt, Y0n(:,3), 'b-', 'LineWidth', 0.9);
ylabel('Rdx'); title('Normal buffering (k=0.10): stable, regular oscillation');
xlim([0 300]); ylim([0 6]); box on;
subplot(2,1,2);
plot(tt, Y0r(:,3), '-', 'Color', [0.85 0.45 0], 'LineWidth', 0.9);
ylabel('Rdx'); xlabel('time');
title('Reduced buffering (k=0.01): chaotic, irregular oscillation');
xlim([0 300]); ylim([0 6]); box on;
print(fig3, 'mito_buffering_timeseries.png', '-dpng', '-r150');
fprintf('Figure saved to mito_buffering_timeseries.png\n');

% ---------------------------------------------------------------
% Local function (must be the last thing in the file: MATLAB requires
% every local function in a script to be defined after all script code)
% ---------------------------------------------------------------
function dx = buffered_rhs(~, x, gA, gB, aA, aB, c, gBuf, aBuf, k, Hact, Hrep)
    % k is the COMBINED Sod1+Sod2 buffering strength for this condition
    Sox=x(1); Prx=x(2); Rdx=x(3); Psi=x(4); Atp=x(5); Suc=x(6); S1=x(7); S2=x(8);
    dx = zeros(8,1);
    dx(1) = -gA*Sox + aA*Hrep(Rdx);
    dx(2) = -gA*Prx + aA*Hrep(Sox);
    dx(3) = -gA*Rdx + aA*Hrep(Prx) + c*Hrep(Suc) - k*(S1+S2)/2*Rdx;
    dx(4) = -gB*Psi + aB*Hrep(Suc);
    dx(5) = -gB*Atp + aB*Hrep(Psi);
    dx(6) = -gB*Suc + aB*Hrep(Atp) + c*Hrep(Rdx) - k*(S1+S2)/2*Suc;
    dx(7) = -gBuf*S1 + aBuf*Hact(Rdx);
    dx(8) = -gBuf*S2 + aBuf*Hact(Rdx);
end

end
