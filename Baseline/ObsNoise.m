function R = ObsNoise(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Process Observation Covariance
% Input:
%       n - Total number of ships
% Ouput:
%       R - Noise covariance (unit: km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = eye(3);
R = blkdiag(I,(1e-5/3600)^2*I,(1e-4)*eye(n-1));