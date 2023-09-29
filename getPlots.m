function figures = getPlots(n,s0,E,a)
% This function calls RK_solver and getThetaNew to plot the phase space,
% parametric, probability density, and phase portrait figures.
% Inputs:
%   n - winding number
%   s0 - initial value of s (usually 0)
%   E - arbitrary guess for energy value
%   a - gamma variable
% Outputs:
%   figures - struct that saves all curves

figures = struct();
[E,s,theta, n, figures.pportait] = RK_solver(n,s0,E,a);

[s,theta] = getThetaNew(s,theta,E,n);

figure
figures.phase = plot(s, theta);
xlabel('s')
ylabel('theta')
grid
xlim([s(1) s(end)])
title({'Saddle Connector for',['n = ', num2str(n)]})

% From the paper, this step relates the probability density to the s and
% theta values
Q = sin(theta);
R_sq = exp(2*cumtrapz(s,Q));
R_sq = R_sq / trapz(R_sq);

figure

figures.prob = plot(s,R_sq);

xlim([s(1) s(end)])
xlabel('s')
ylabel('Probability Density R^2')
title({'Probability density plots for' ,['n = ', num2str(n)]}')
xlim([-20 20])  % creating plotting x-limits to capture characteristics


% Utilizing the Prufer Transform to relate the probability density to the
% parametric variables
R = sqrt(R_sq);
u = R.*cos(theta/2);
v = R.*sin(theta/2);

figure

figures.parametric = plot(u,v);
xlabel('u')
ylabel('v')
title({'Parametric plot in u-v space:', ['n = ', num2str(n)]})

end
