function [mu,sigma] = invUT(Xsigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inversed Unscented Transformation
%Input: Sigma point - Vector
%Output: (mu, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 2;
n = (size(Xsigma,2)-1)/2; %Size of state = M*N

w0 = lambda/(lambda+n);
wi = 1/(2*(lambda+n));
mu = w0*Xsigma(:,1);
for i = 1:2*n
    mu = mu+wi*Xsigma(:,i+1);
end

sigma = w0*(Xsigma(:,1)-mu)*(Xsigma(:,1)-mu).';
for i = 1:2*n
    sigma = sigma+wi*(Xsigma(:,i+1)-mu)*(Xsigma(:,i+1)-mu).';
end