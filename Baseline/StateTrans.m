function X1 = StateTrans(X0,delT,T,LunarData)
%%%%%%%%%%%%%%%% 
% Compute ECI Dynamics with RK4 (Units in km or km/s)
% Input: 
%       X0 - ECI State Vector at Current Timestep: [r0, v0, r1, v1,..., rn, vn].'
%       delT - Time Step
%       T - Current Real Time
%       LunarData - Lunar Data Table
%Output:
%       X1 - ECI State Vector at Next Timestep: [r0, v0, r1, v1,..., rn, vn].'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Interpolated data at t, t+h/2, t+h (Same units as Lunar Data: km & km/s)
MVec1 = LunarDataInterp(T,LunarData);
MVec2 = LunarDataInterp(T+delT/2,LunarData);
MVec3 = LunarDataInterp(T+delT,LunarData);

delTm = seconds(delT); %Convert any timestep value to seconds
k1 = delTm*Dynamics(X0,MVec1);
k2 = delTm*Dynamics(X0+k1/2,MVec2);
k3 = delTm*Dynamics(X0+k2/2,MVec2);
k4 = delTm*Dynamics(X0+k3,MVec3);

X1 = X0 + (k1+2*k2+3*k3+k4)/6; %Runge-Kutta 4th Order Method