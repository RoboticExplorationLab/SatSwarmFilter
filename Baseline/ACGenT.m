function [A,C] = ACGenT(T)
%%%%%%%%%%%%%%%% 
% Generate Linearized A, C at Time T
% Input: 
%       T - Current Real Time
% Output:
%       A, C Matrix at Time T
% Comment:
%       Here we used input T directly for the computation of A&C.
%       Meanwhile, we also used the true data as the staring point for us
%       to apply RK4 to get the current ship and lunar positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = datetime(T,'ConvertFrom','datenum'); %Convert to datetime format
[x,MVec] = PosGen(T); %Compute current vector position 
[A,C] = ACGenPos(x,MVec); %Compute A & C matrices