function [x,MVec] = PosGen(T)
%%%%%%%%%%%%%%%%
% Generate ECI position for both the swarm and the lunar at time T
% Input: 
%       T - Current Real Time
% Output:
%       x - ECI Positions for swarm
%       MVec - ECI Position for lunar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute current vector position 
% T = datetime(T,'ConvertFrom','datenum'); %Convert to datetime format
global Xtrue T_SIM LunarData

idxI = find(T>=T_SIM,1,'last'); %Index of last sampled time 

T0 = T_SIM(idxI); %Last sampled time
X0 = Xtrue(:,idxI); %Position at last sampled time
if isempty(T-T0)
    x = X0; %If current time is sampled time
else
    x = StateTrans(X0,T-T0,T0,LunarData); %ECI Vector at Time T
end
MVec = LunarDataInterp(T,LunarData); %Lunar ECI Vector at Time T