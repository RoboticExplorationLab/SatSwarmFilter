function Q = ProcessNoise(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Process Noise Covariance
% Input:
%       n - Total number of ships
% Ouput:
%       Q - Noise covariance (unit: km or km/s)
% Comment:
%       The state uncertainty is 1km for position and 1cm/s for velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qi = [ones(1,3),(1e-4)^2*ones(1,3)];
Q = diag(repmat(Qi,1,n));