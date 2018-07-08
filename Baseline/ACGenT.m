function [A,C] = ACGenT(T)
%%%%%%%%%%%%%%%% 
% Generate Linearized A, C at Time T
% Input: 
%       T - Current Real Time in datenum format
%       Xtrue - ECI Positions at Time Points (treated as constants)
%       T_SIM - Time Vecotr correspondsing to Xtrue (treated as constants)
% Output:
%       A, C Matrix at Time T
% Comment:
%       Here we used input T directly for the computation of A&C.
%       Meanwhile, we also used the true data as the staring point for us
%       to apply RK4 to get the current ship and lunar positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute current vector position 
% T = datetime(T,'ConvertFrom','datenum'); %Convert to datetime format
global Xtrue T_SIM LunarData

idxI = find(T>=T_SIM,1,'last'); %Index of last sampled time 

T0 = T_SIM(idxI); %Last sampled time
X0 = Xtrue(:,idxI); %Position at last sampled time
if isempty(T-T0)
    X1 = X0; %If current time is sampled time
else
    X1 = StateTrans(X0,T-T0,T0,LunarData); %ECI Vector at Time T
end
MVec = LunarDataInterp(T,LunarData); %Lunar ECI Vector at Time T

[A,C] = ACGenPos(X1,MVec);