function [ w ] = statdist( P )
%STATDISTBIO Returns stationary distribution of biotic states
%   For use with Shannon-Wiener index calculation
%
%INPUT
% P - transition matrix
% n - index for bare space

[l,W,V] = sorted_eigs(P);

%stationary distribution 
w = W(:,1);



end

