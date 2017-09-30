function [ dwn, dwnscale ] = Pstat_partialderiv( P,i,j,n )
%PSTAT_PARTIALDERIV calculates the partial derivative of the nth eigenvector
%   Partial w_n over partial p_ij (or i,j element of the P matrix)
% SOURCE
%       Matthew Spencer (2006), https://www.liverpool.ac.uk/~matts/markov.html
%       with slight modifications
% INPUT
%   P - column stochastic matrix; element p_ij is the probability of i given j
%   i,j - index of the P element being perturbed
%   n - index of the eigenvector whose sensitivity we want to evaluate
% OUTPUT
%   dwn - vector; partial w_n over partial p_ij
%   dwnscale - partial w_n / ||w_n|| over partial p_ij

s = size(P,1);
[lambda,W,V] = sorted_eigs(P);

%get the complex conjugates of the left eigenvectors
Vcc = conj(V);   

m = [1:s];
m(n) = [];   %take out the nth index

%do the summation
dwn = 0;
for k = 1:(s-1)
    %Spencer
    %dwn=dwn+((Vcc(m(k),i)-Vcc(m(k),j))/(lambda(n)-lambda(m(k))))*W(:,m(k));
    %Hill
    dwn = dwn + (Vcc(m(k),i) / (lambda(n) - lambda(m(k)))) * W(:,m(k));
end

dwn = W(j,n) * dwn;

%get scaled version
dwnscale = dwn - W(:,n)*sum(dwn);

end

