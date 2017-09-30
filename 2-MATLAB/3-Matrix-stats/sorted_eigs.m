function [ lambda,W,V ] = sorted_eigs( P )
%SORTED_EIGS returns eigenvalues and eigenvectors
%   sorted according to descending real part of the eigenvalue
% SOURCE
%   Matthew Spencer (2006), https://www.liverpool.ac.uk/~matts/markov.html
% INPUT
%   P - a discrete-time matrix, P
% OUTPUT
%   lambda - a vector of eigenvalues, 
%          - sorted in descending order of real part magnitude
%   W - matrix of corresponding right eigenvectors (columns)
%     - w1 (first column) scaled to sum to 1
%   V - matrix of corresponding left eigenvectors (rows)

% calc. eigenvalues and eigenvectors
[W,d] = eig(P);  %d is a diagonal mx of lambdas
lambda = diag(d);

% sort eigenvalue and eigenvectors according to lambda dominance (real part only)
[ldum,ind] = sort(real(lambda),'descend');
lambda = lambda(ind);
W = W(:,ind);
W = W/sum(W(:,1));  %scale so that dominant eigenvector sums to 1
V = conj(inv(W));   %rows of V are left eigenvectors

end

