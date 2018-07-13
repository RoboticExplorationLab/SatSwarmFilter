function [f1, C1] = NumIntA(Tmin, Tmax, At)
%%%%%%%%%%%%%%%% 
% Numerically integrate A and returns C matrix at Tmin
% Input: 
%       A - A(t) function takes datetime input and output A matrix
%       tmin, tmax - Start and ending of integration
% Output:
%       f - Integral of A(t) from tmin to tmax
% Comment:
%       Here we used input T directly for the computation of A&C.
%       Meanwhile, we also used the true data as the staring point for us
%       to apply RK4 to get the current ship and lunar positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 100;
Tvec = linspace(Tmin, Tmax, n); %Range vector (datetime format)
Tstep = (Tvec(2) - Tvec(1)); %Time step in seconds (datetime format)
TstepNum = seconds(Tstep); %Time step in seconds (num format)
[A1, C1] = At(Tmax);
f1 = zeros(size(A1));
for i = 1:n-1
    T = Tvec(i)+ Tstep/2; %Midpoint quadrature rule (datetime format)
    [A,~] = At(T); %Compute A matrix (1/s)
    f1 = f1 + A*TstepNum; %Integration (num format)
end