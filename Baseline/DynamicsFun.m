function xdot = DynamicsFun(x,MVec)
%%%%%%%%%%%%%%%% 
% Derivative of ECI Vector
% Input: 
%       x - ECI State Vector at Current Timestep: [r0, v0, r1, v1,..., rn, vn].'
%       Mvec - ECI Vector for Luna at Current Time
% Output:
%       xdot - Derivative of ECI State Vector at Current Timestep
% Comment:
%       Here we used input t implicitly in Mvec since it is solely
%       determined by time. This vector should be computed with time in the
%       mother function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global muE muM %km^3/s^2
if isrow(MVec)
    MVec = MVec.';
end
m = 3; %Coordinate Size
n = length(x)/(2*m); %Number of targets

for idx = 1:n %Isolate each target separately
    rivec = x(2*m*(idx-1)+1 : 2*m*(idx-1)+m); %Position vector in km
    vivec = x(2*m*(idx-1)+m+1 : 2*m*idx,:); %Velocity vector in km/s
    xdot(2*m*(idx-1)+1 : 2*m*(idx-1)+m,:) = vivec; %Position derivative in km/s
    
    rMvec = rivec - MVec(1:3); %Vector from moon to satellite
    rMnorm = norm(rMvec); %Distance between moon and satellite
    rEnorm = norm(rivec); %Distance between earth and satellite
    
    avec = -muM/(rMnorm^3)*rMvec - muE/(rEnorm^3)*rivec; %Newton's 2nd Law
    
    xdot(2*m*(idx-1)+m+1 : 2*m*idx,:) = avec; %Velocity derivative
end