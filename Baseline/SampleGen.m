function X0 = SampleGen(OE0,M)
%%%%%%%%%%%%%%%% 
% Generate sample ships
% Input: 
%       OE0 - Initial OE of Hub
%       M - Total Number of Ships
% Output:
%       X0 - Initial ECI Vector of Children (1st cell open for hub)
% Comment:
%       Here we vary the semi-major axis by generating a uniform sample
%       from [.5, .8] and then multiplying by 9000km. The final initial
%       variation between different children and hub is less than
%       1000 km WHP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X0{1} = s3_state_keptocart(OE0, 1)/1e3;
M = M-1; %Number of ships excepting huvb
for i = 1:M
    OE = OE0;
    OE(1)=OE(1) +(rand*.3+.5)*9000e3; %Vary semi-major axis randomly (in m)
    X0{i+1} = s3_state_keptocart(OE, 1)/1e3; %Initial ECI vector of children i in km or km/s
end