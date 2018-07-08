function f2 = NumIntA1(Tmin,Tmax,At)
%%%%%%%%%%%%%%%% 
% Numerically integrate the second order term of Peano-Baker series
% \int_{t_{1}}^{t_{2}}A(\tau_{1})\int_{t_{1}}^{\tao_{1}}A(\tao_{2})d\tao_{1}d\tao{2}
% Input: 
%       A - A(t) function takes datetime input and output A matrix
%       tmin, tmax - Start and ending of integration
% Output:
%       f1 - Integral
% Comment:
%       This is a doouble integral, thus we will use NumIntA to integrate A
%       and this function for the entire integral. The inner integral uses
%       100 points (NumIntA). Outer integral uses 10 points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 10; %Number of outer integral points
Tvec = linspace(Tmin, Tmax, n); %Range vector (datetime format) for outer integral
Tstep = (Tvec(2) - Tvec(1)); %Time step in seconds (datetime format) for outer integral
TstepNum = seconds(Tstep); %Time step in seconds (num format) for outer integral
[A0, ~] = At(Tmin);
f2 = zeros(size(A0*A0));
for i = 1:n-1
    Aint = NumIntA(Tmin, Tvec(i), At); %Inner integral
    T = Tvec(i)+ Tstep/2; %Midpoint quadrature rule (datetime format)
    [A,~] = At(T); %Compute A matrix (1/s)
    f2 = f2 + A*Aint*TstepNum; %Integration (num format)
end