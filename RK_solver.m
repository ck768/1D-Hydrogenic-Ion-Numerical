function [E,s,theta,n, pportrait] = RK_solver(n,s0,E,a)
% This function solves the nonlinear ODE system using the Runge-Kutta 4
% method. It implements a binary search that is bounded by E=1 and E=-1,
% and then iterates through the binary search to find the saddle connector
% solution.
% Inputs:
%   n - winding number
%   s0 - initial value of s (usually 0)
%   E - arbitrary guess for energy value
%   a - gamma variable
% Outputs:
%   E - energy value for saddle connector solution
%   s - s-array
%   theta - theta-array
%   n - winding number
%   pportrait - phase portrait plot


% We have proven in our paper that for s=0, the initial value of Theta is
% equal to -n*pi. We have also proven symmetry in our system, so we solve
% only the system for s=0 to s=inf, and then reflect about the y-axis.
theta0 = -n*pi;

Eb = 1;
Ea = -1;
tol = 1e-15;

figure
while abs(Ea-Eb)>tol
    
    h = 0.01; % time step
    N = 10000; % number of time steps
    
    s = zeros(1, N+1); % array to store s values
    theta = zeros(1, N+1); % array to store theta values
    s(1) = s0; % set initial s value
    theta(1) = theta0; % set initial theta value
    
    for i = 1:N
        k1 = h*diff_eq(s(i), theta(i), a, E);
        k2 = h*diff_eq(s(i)+h/2, theta(i)+k1/2, a, E);
        k3 = h*diff_eq(s(i)+h/2, theta(i)+k2/2, a, E);
        k4 = h*diff_eq(s(i)+h, theta(i)+k3, a, E);
        theta(i+1) = theta(i) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        s(i+1) = s(i) + h;
    end
    
    theta_final = theta(end);
    
        if theta_final < -2*n*pi - acos(E)                                                                                                                 % the  lower half of the E-interval., otherwise take the upper half
            Eb = E;
        else
            Ea = E;
        end
        E = (Ea+Eb)/2; 
    plot(s,theta);
    hold on
end
pportrait = plot(s,theta);
xlim([s(1) s(end)])
xlabel('s')
ylabel('\Theta')
title('\Theta - s Phase Portrait')

%% local function - our nonlinear ODE system
function dtheta = diff_eq(s, theta, a, E)
    dtheta = 2*cos(theta) - a*exp(-abs(s)) - 2*E;
end
end