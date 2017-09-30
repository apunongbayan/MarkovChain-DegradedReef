function [ dPmats ] = sensitivities( Psamples, weight )
%SENSITIVITIES Returns sensitivity matrices calculated for each P sample
%
% Input
%   Psamples - transition probabilities, sample = row vector
%   weight - how to weight the state vector; default: only hard corals
% Output
%   dPmats - 3 dim. matrix; z - samples; x,y - {sensitivity of w to pij}

if nargin < 2
    weight = [0 0 0 0 0 1 1 0];
end
     
compens = 3;                        % compensation type - proportional
[nsamps,nparms] = size(Psamples);   
s = sqrt(nparms);                   % number of states

dPmats = zeros(s,s,nsamps);
for ni = 1:nsamps
    Phat = reshape(Psamples(ni,:),8,8);
    dPmat = zeros(s,s);
    for i = 1:s
        for j = 1:s
            [sens,elas] = sens_statdistcarbonate(Phat,i,j,weight,compens);
            dPmat(i,j) = dPmat(i,j) + sens;
        end
    end
    %dPmat
    %size(dPmat)
    dPmats(1:s,1:s,ni) = dPmat;
    ni
end

end

