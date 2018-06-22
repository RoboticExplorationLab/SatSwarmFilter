function Q = ProcessNoise(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Process Noise Covariance
% Input:
%       n - Total number of ships
% Ouput:
%       Q - Noise covariance (unit: km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 1e-12*eye(6*n);