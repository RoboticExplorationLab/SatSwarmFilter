function R = ObsNoise(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Process Observation Covariance
% Input:
%       n - Total number of ships
% Ouput:
%       R - Noise covariance (unit: km or km/s)
% Comment:
%       The observation uncertainty is 1km for hub position, 1cm/hr for
%       range rate, 10 m for relative range.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = eye(3);
R = blkdiag(I,(1e-5/3600)^2*I,(1e-4)*eye(n-1));