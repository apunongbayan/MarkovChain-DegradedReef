function [ allProps, propNames ] = communityProperties( jagssamples, ...
   resoln, eindex )
%COMMUNITYPROPERTIES Compiles successional indices' values calculated from
%   random samples of the transition probability matrix
%   depends on calc_communityProps.m
%
% Input
%   jagssamples - a matrix; row - sample, col - transition probability
%   resoln - of the transition probabilities; in months
%   eindex - index of empty state (CCA + turf)
% Output
%   allProps - a cell array containing
%       convergence measures
%       recurrence times
%       turnover times
%       normalized entropy
%       etc. (as descibed in Hill et al. 2004)
%   propNames - the name of the variable in each cell

if nargin < 3
    resoln = 1;
    eindex = 1;
end

dampvec = [];
rectimesvec = [];  %rows - states; cols - samples; in years
turntimesvec = [];
brrectimevec = []; %in years
normentrovec = [];
dominantvec = [];
dobcoeffvec = [];       % 0-converges to stat dist in 1 iter; 1-no converg.
convergeratevec = [];   %x100 percent per month'
halflifevec = [];       %months
nlogdobcoeffvec = [];   %x100 percent per month
repby_jvec = [];
entropycolsvec = [];    %entropy in functional group transitions

[ns, np] = size(jagssamples);

for i = 1:ns
    %phat = reshape(jagssamples(randi(ns,1),:),8,8)
    phat = reshape(jagssamples(i,:),8,8);  
    [ repby_jvec(:,end+1), normentrovec(end+1), turntimesvec(:,end+1),...
        rectimesvec(:,end+1), brrectimevec(end+1),...
        entropycolsvec(:,end+1), dobcoeffvec(end+1),...
        dampvec(end+1), convergeratevec(end+1), halflifevec(end+1),...
        nlogdobcoeffvec(end+1),dominantvec(:, end+1)] ...
        = calc_communityProps( phat, eindex, resoln);   
end

propNames = {'Damping ratio','Recurrence times (years)',...
    'Turnover times (years)','Clear space recurrence time',...
    'Normalized entropy','Stationary distribution',...
    'Dobrushin coefficient','Convergence rate','Half-life',...
    'n log (Dobrushin coeff)','repby_jvec',...
    'Entropy of transitions given state'};
allProps{1} = dampvec;          
allProps{2} = rectimesvec;
allProps{3} = turntimesvec;
allProps{4} = brrectimevec;
allProps{5} = normentrovec;
allProps{6} = dominantvec;
allProps{7} = dobcoeffvec;
allProps{8} = convergeratevec;
allProps{9} = halflifevec;
allProps{10} = nlogdobcoeffvec;
allProps{11} = repby_jvec;
allProps{12} = entropycolsvec;

end

