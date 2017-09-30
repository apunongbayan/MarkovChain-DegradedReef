function [ Pdwn, Pdwnscale ] = Pstat_deriv( P, k, V )
%PSTAT_DERIV returns the total derivative of the kth eigenvector
%   of the P matrix wrt to a compensation pattern V
% SOURCE
%   Matthew Spencer (2006), https://www.liverpool.ac.uk/~matts/markov.html
% INPUT
%   P - matrix of transition probabilities
%   k - index of the eigenvector to evaluate
%   V - matrix of compensation patterns
% OUTPUT
%   Pdwn - col vec; total derivative of the nth eigenvector wrt to perturbations in P given by V
%   Pdwnscale - total derivative of the scaled nth eigenvector...

% Note V is a matrix of zeros, and non-zero entries depend on the compensation pattern
% If SPECIFIC:
%    dp_ij = 1; dp_jj = -1
% If UNIFORM
%    dp_ij = 1; dp_mj = - 1 / (s-1)
% If PROPORTIONAL
%    dp_ij = 1; dp_mj = - p_mj / (1 - p_ij)


s = size(P,1);
Pdwn = zeros(s,1);
Pdwnscale = Pdwn;

for m=1:s
    for n=1:s
        [dwn,dwnscale] = Pstat_partialderiv(P,m,n,k); %kth evec
        Pdwn = Pdwn + V(m,n)*dwn;   
        Pdwnscale = Pdwnscale + V(m,n)*dwnscale;
    end
end

