% ADAPTIVE_LORENZ_DEMO
%
% A 4-dimensional nonlinear system (x, y, z, k) built from the classic
% 3D Lorenz equations plus one extra "regulatory" dimension, k, that
% adjusts a feedback gain.
%
%   dx/dt = sigma*(y - x)
%   dy/dt = x*(rho - z) - y - k*y
%   dz/dt = x*y - beta*z
%   dk/dt = gamma*y^2
%
% sigma, rho, beta are fixed at the classic chaotic Lorenz values
% (10, 28, 8/3), so the bare (x,y,z) subsystem on its own is the
% textbook chaotic attractor with a positive Lyapunov exponent and
% the familiar sensitivity to initial conditions.
%
% k has no dynamics of its own beyond integrating y^2, so it never
% decreases, and it only stops changing once y is driven toward zero.
% It behaves like a homeostat: any sustained oscillation in y makes k
% climb, feeding more damping back into the y-equation, until the
% oscillation is suppressed.
%
% Two regimes:
%   FREE     - all four dimensions (x,y,z,k) co-evolve. k climbs and
%              drives the trajectory to the origin, an equilibrium
%              that always exists for the Lorenz equations. Two
%              trajectories started a tiny distance apart converge to
%              indistinguishable, regardless of where they started:
%              robust to initial conditions.
%   CLAMPED  - k is held fixed at its initial value instead of being
%              allowed to adapt, the analogue of clamping a regulatory
%              variable in a reduced/in-vitro preparation instead of
%              letting it respond to the rest of the system. With k
%              pinned at 0, the remaining (x,y,z) subsystem IS the
%              standard chaotic Lorenz attractor: nearby trajectories
%              separate exponentially.
%
% A scan over fixed values of k shows the transition is sharp: any
% appreciable nonzero fixed damping already restabilizes the system,
% so it specifically matters that the clamp catches the regulatory
% variable at its inactive/off value, not merely that it is held
% fixed at all.
%
% Run this file directly (MATLAB or GNU Octave) to reproduce the
% printed numbers and the divergence plot (saved as divergence.png).

function adaptive_lorenz_demo()

    sigma = 10.0;
    rho   = 28.0;
    beta  = 8.0/3.0;
    gamma = 8.0;   % adaptation rate for the regulatory dimension k

    eps  = 1e-8;
    x0   = 1.0;
    tEnd = 40.0;
    nPts = 6000;
    tSpan = linspace(0, tEnd, nPts);

    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

    % --- FREE system: (x,y,z,k) all evolve together ---
    freeRHS = @(t, s) [ sigma*(s(2) - s(1)); ...
                        s(1)*(rho - s(3)) - s(2) - s(4)*s(2); ...
                        s(1)*s(2) - beta*s(3); ...
                        gamma*s(2)^2 ];

    [~, Sa] = ode45(freeRHS, tSpan, [x0,       1.0, 1.0, 0.0], opts);
    [~, Sb] = ode45(freeRHS, tSpan, [x0 + eps, 1.0, 1.0, 0.0], opts);
    distFree = sqrt(sum((Sa(:,1:3) - Sb(:,1:3)).^2, 2));

    % --- CLAMPED system: k fixed at 0, only (x,y,z) evolve ---
    kFixed = 0.0;
    clampedRHS = @(t, s, k) [ sigma*(s(2) - s(1)); ...
                              s(1)*(rho - s(3)) - s(2) - k*s(2); ...
                              s(1)*s(2) - beta*s(3) ];

    [~, Sc] = ode45(@(t,s) clampedRHS(t, s, kFixed), tSpan, [x0,       1.0, 1.0], opts);
    [~, Sd] = ode45(@(t,s) clampedRHS(t, s, kFixed), tSpan, [x0 + eps, 1.0, 1.0], opts);
    distClamped = sqrt(sum((Sc - Sd).^2, 2));

    fprintf('FREE system   : separation grows from %.1e to %.1e (k settles at %.1f) -> converges, robust to initial conditions\n', ...
        distFree(1), distFree(end), Sa(end,4));
    fprintf('CLAMPED (k=0) : separation grows from %.1e to %.1e -> diverges, sensitive to initial conditions (chaos)\n', ...
        distClamped(1), distClamped(end));

    fprintf('\nScan of fixed k values (clamped system) -- how sharp is the transition?\n');
    kScan = [0, 0.5, 1, 2, 5, 10, 20];
    for kk = kScan
        [~, Se] = ode45(@(t,s) clampedRHS(t, s, kk), tSpan, [x0,       1.0, 1.0], opts);
        [~, Sf] = ode45(@(t,s) clampedRHS(t, s, kk), tSpan, [x0 + eps, 1.0, 1.0], opts);
        dScan = sqrt(sum((Se - Sf).^2, 2));
        fprintf('  k=%5.1f  final separation = %.2e\n', kk, dScan(end));
    end

    % --- Plot ---
    fig = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');

    subplot(1, 2, 1);
    semilogy(tSpan, distFree, 'Color', [0.169, 0.424, 0.690], 'LineWidth', 1.8); hold on;
    semilogy(tSpan, distClamped, 'Color', [0.753, 0.224, 0.169], 'LineWidth', 1.8);
    xlabel('time');
    ylabel('distance between two nearby trajectories');
    title('Separation of initially close trajectories (initial offset = 1e-8)');
    legend('all dimensions free (x,y,z,k)', 'k clamped at 0 (x,y,z only)', ...
           'Location', 'east', 'FontSize', 8);
    ylim([1e-19, 1e2]);
    grid on;

    subplot(1, 2, 2);
    plot(Sc(:,1), Sc(:,3), 'Color', [0.753, 0.224, 0.169], 'LineWidth', 0.5);
    title('Clamped system (k=0): bare chaotic Lorenz attractor');
    xlabel('x');
    ylabel('z');

    saveas(fig, 'divergence.png');
    fprintf('\nSaved divergence.png\n');

end
