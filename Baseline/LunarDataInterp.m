function MVec = LunarDataInterp(t,LunarData)
%%%%%%%%%%%%%%%% 
% Linear Interpolate Lunar Data at Current Time
% Input: 
%       t - Current Real Time
%       LunarData - Lunar Data Table
%Output:
%       MVec - Interpolated Value of ECI Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxF = find(t<LunarData{1},1); %Index of next lunar data time 
idxI = find(t>LunarData{1},1,'last'); %Index of last lunar data time 

MTime = [LunarData{1}(idxI),LunarData{1}(idxF)].'; %Lunar data time vector
MVec0 = [LunarData{2}(idxI,:);LunarData{2}(idxF,:)]; %Lunar data vector

MVec = interp1(MTime,MVec0,t); 