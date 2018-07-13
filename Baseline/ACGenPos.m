function [A,C] = ACGenPos(x,MVec)
%%%%%%%%%%%%%%%% 
% Generated Linearized A, C matrix using vector input
% Input: 
%       x - ECI State Vector at Current Timestep: [r0, v0, r1, v1,..., rn, vn].'
%       Mvec - ECI Vector for Luna at Current Time
% Output:
%       A, C Matrix
% Comment:
%       Here we used input t implicitly in Mvec and x since it is solely
%       determined by time. This vector should be computed with time in the
%       mother function. It uses positions as inputs with time being
%       embedded implicitly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global muE muM %km^3/s^2
if isrow(MVec)
    MVec = MVec.';
end
m = 3;
n = length(x)/(2*m); %Number of targets
I = eye(m); %Identity matrix
Z = zeros(m); %Zero matrix fo A
r0 = x(1:m); %Hub position

TotL = length(x); %Total width of the C matrix for zero padding (6 for identity matrix)
C = blkdiag(eye(3),eye(3)); %Hub observation derivative
C = padarray(C,[0,TotL - size(C,2)],'post');

for idx = 1:n
    rivec = x(2*m*(idx-1)+1 : 2*m*(idx-1)+m); %Position vector in km
    vivec = x(2*m*(idx-1)+m+1 : 2*m*idx,:); %Velocity vector in km/s
    
    rMvec = rivec - MVec(1:3); %Vector from moon to satellite
    rMnorm = norm(rMvec); %Distance between moon and satellite
    rEnorm = norm(rivec); %Distance between earth and satellite
    
    %A Computation
    dA1 = muM*(3*rivec*(rivec-rMvec).'-I*rMnorm^2)/rMnorm^5;
    dA2 = -3*muM*rMvec*(rivec-rMvec).'/rMnorm^5;
    dA3 = muE*I/rEnorm^3 - 3*muE*(rivec*rivec.')/rEnorm^5;
    dAi = dA1 +dA2 + dA3; %dAi
    dAi = dAi*3600; %Convert to min^-2
    blkAi{idx} = [Z,I;dAi,Z]; %Block matrices
    
    %C computation
    if C ~= 1
        dCi1 = (r9 - rivec).'/norm(rivec - r0); %d/dr0||ri-r0||
        dCi2 = (rivec - r0).'/norm(rivec - r0); %d/dri||ri-r0||
        dCiRow = [dCi1,zeros(1,6*idx-9),dCi2];
        CiRow = padarray(dCiRow,[0,TotL - size(dCiRow,2)],'post');
        C = [C;CiRow];
    end
end
A = blkdiag(blkAi{:});