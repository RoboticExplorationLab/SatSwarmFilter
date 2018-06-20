function X1 = StateTrans(X0,delT)
%%%%%%%%%%%%%%%% Compute Dynamics of ROE
%Input: 
%       X - ECI State Vector at Current Timestep: [r0, v0, r1, v1,..., rn, vn].'
%       Q - Process Noise Covariance
%Output:
%       X1 - ECI State Vector at Next Timestep: [r0, v0, r1, v1,..., rn, vn].'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = s3_constants;
mu = const.GM.EARTH;

m = 3; %Coordinate Size
n = length(X0)/(2*m); %Number of targets

for idx = 1:n
    r0vec = X0(2*m*(idx-1)+1 : 2*m*(idx-1)+m); %Position vector
    v0vec = X0(2*m*(idx-1)+m+1 : 2*m*idx); %Velocity vector
    
    r0norm = norm(r0vec);
    avec = mu/(r0norm^3)*r0vec;
    v1vec = v0vec + avec*delT;
    r1vec = r0vec + v0vec*delT + avec/2*delT^2;
    
    X1(2*m*(idx-1)+m+1 : 2*m*idx) = v1vec; %Update velocity
    X1(2*m*(idx-1)+1 : 2*m*(idx-1)+m) = r1vec; %Update position
end