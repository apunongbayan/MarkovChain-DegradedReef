function [ x0 ] = sampleInitDistfromAlpha( alpha0 )
%SAMPLEINITDISTFROMALPHA Sample initial probability vector from Dirichlet
%prior
% Input
%   alpha0 - a vector; 
% Output
%   x0 - vector for initial state distribution

theta_raw = gamrnd(alpha0,1);
x0 = theta_raw / sum(theta_raw);


end
