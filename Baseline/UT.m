function [xi, wi] = UT(mu, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unscented Transformation
%Input: (mu, sigma)
%Output: Sigma point (xi, wi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 2; %User defined lambda 
n = length(mu); %Size of state = M*N

xi(:,1) = mu;
wi(1) = lambda/(n+lambda);
wi(2:2*n+1) =1/(2*(n+lambda));

Mat = sqrtm((n+lambda)*sigma);
for i = 1:n
    xi(:,i+1) = mu + Mat(:,i);    
end
for i = n+1:2*n
    xi(:,i+1) = mu - Mat(:,i-n);
end
