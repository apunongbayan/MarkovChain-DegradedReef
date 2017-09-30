function [ Psample ] = samplePmatfromAlpha( Alpha )
%SAMPLEPMATFROMALPHA Sample a probability matrix from Dirichlet matrix

K = size(Alpha,1);

%Build new P matrix
Pp = zeros(K,K);

for j = 1:K
    theta_raw = gamrnd(Alpha(:,j),1);
    theta = theta_raw / sum(theta_raw);
    Pp(:,j) = theta;
end

Psample = Pp;


end

