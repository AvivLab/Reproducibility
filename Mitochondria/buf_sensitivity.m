function buf_sensitivity()
% Parameter scan supporting mito_redox_buffering.m: single-buffer version of
% the same two-ring cross-talk network (Os-Sk-Hg / Ca2-Cn-Ri, cross-talk
% Hg<->Ri), used only to map out how a single combined buffering-strength
% parameter k affects the twin-trajectory separation at t=250, across
% k = 0 to 0.15. Saves the scan to sensitivity_data.mat, which
% mito_redox_buffering.m loads (if present) to draw the sensitivity figure.

n=4; K=1;
Hact = @(x) (max(x,0).^n) ./ (K^n + max(x,0).^n);
Hrep = @(x) 1 - Hact(x);
gA=1.0; gB=2.0; aA=8.0; aB=8.0; c=3.0; gBuf=1.5; aBuf=8.0;

opts = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',1.0);
Tspan=[0 250]; tt=linspace(0,250,2000);
x0 = [0.42 0.72 0.10 0.68 0.24 0.19 0.30]'; x1=x0; x1(1)=x1(1)+1e-6;

kvals = 0:0.01:0.15;
seps = zeros(size(kvals));
for i = 1:length(kvals)
  k = kvals(i);
  rhs = @(t,x) single_buf_rhs(t,x,gA,gB,aA,aB,c,gBuf,aBuf,k,Hact,Hrep);
  [t0,y0]=ode45(rhs,Tspan,x0,opts);
  [t1,y1]=ode45(rhs,Tspan,x1,opts);
  Y0=interp1(t0,y0,tt); Y1=interp1(t1,y1,tt);
  sep=sqrt(sum((Y0-Y1).^2,2));
  seps(i) = sep(end);
  fprintf('k=%.3f | sep(end)=%.6f\n', k, sep(end));
end
save('sensitivity_data.mat', 'kvals', 'seps');

end

% ---------------------------------------------------------------
% Local function (must be the last thing in the file: MATLAB requires
% every local function in a script to be defined after all script code;
% keeping the same rule here even though this file is now a function file,
% where it is not strictly required, for consistency with the other
% scripts in this project)
% ---------------------------------------------------------------
function dx = single_buf_rhs(~, x, gA, gB, aA, aB, c, gBuf, aBuf, k, Hact, Hrep)
    Os=x(1); Sk=x(2); Hg=x(3); Ca2=x(4); Cn=x(5); Ri=x(6); Bf=x(7);
    dx = zeros(7,1);
    dx(1) = -gA*Os  + aA*Hrep(Hg);
    dx(2) = -gA*Sk  + aA*Hrep(Os);
    dx(3) = -gA*Hg  + aA*Hrep(Sk) + c*Hrep(Ri) - k*Bf*Hg;
    dx(4) = -gB*Ca2 + aB*Hrep(Ri);
    dx(5) = -gB*Cn  + aB*Hrep(Ca2);
    dx(6) = -gB*Ri  + aB*Hrep(Cn) + c*Hrep(Hg) - k*Bf*Ri;
    dx(7) = -gBuf*Bf + aBuf*Hact(Hg);
end
