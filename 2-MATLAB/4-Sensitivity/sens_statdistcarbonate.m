function [ sens, elas ] = sens_statdistcarbonate( P, i, j, c, ctype )
%SENS_STATDISTCARBONATE Returns sensitivity /elasticity of stationary 
%   distribution of carbonate states to pijs
%
%INPUT
% P - transition matrix
% i,j - element of P being perturbed
% c - density weights
% ctype - compensation pattern to maintain column stochasticity of P
%       1 - Specific: dp_ij = 1; dp_jj = -1
%       2 - Uniform: dp_ij = 1; dp_mj = - 1 / (s-1)
%       3 - Proportional
%
%OUTPUT
% sens - sensitivity of weighted distribution to change in pij
% elas - proportional sensitivity of weighted dist'n to prop change in pij

[l,W,V] = sorted_eigs(P);   %[eigenvals, right evecs, left evec
w = W(:,1);                 %take dominant right evec
nstates = size(P,1);

%sensitivity stationary distribution (all states)
dp = getCompensation(P,i,j,ctype);
[dw1_unscaled, dw1] = Pstat_deriv(P,1,dp);  %get total derivative of 1st
% eigenvector of P wrt to compensation pattern dp

%indices of states / how densities should be weighted
%example: c = [0 0 0 0 0 1 1 0] -- only scleractinian states
idum = find(ne(c,0));
sens = sum(dw1(idum));
elas = (P(i,j)/sum(w(idum))) * sens;    %scaled according to stationary 
% density of selected states

end

