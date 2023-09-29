function [s,theta_new] = getThetaNew(s,theta,E,n)
% This function finds where the saddle connector solution is expected to 
% end at (through analytical methods) and then creates a modified 
% theta-array that tracks when the theta value is close to the expected.
% Inputs from RK_solver:
%   s - s-array (from s=0 to s=inf)
%   theta - theta-array
%   E - energy value for saddle connector solution
%   n - winding number
% Outputs:
%   s - s-array (from s=-inf to s=inf)
%   theta_new - theta-array with corrected end behavior

theta_new = theta;
ind = find(abs(theta - (-2*n*pi-acos(E))) < 1e-4,1);
theta_new(ind:end) = -2*n*pi-acos(E);

% This next step mirrors the curve about the y-axis due to symmetry
theta_new = [-(flip(theta_new)+2*n*pi) theta_new];
s = [-flip(s) s];
end
