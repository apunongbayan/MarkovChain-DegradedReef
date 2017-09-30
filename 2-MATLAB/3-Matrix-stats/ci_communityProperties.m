function [ muArray, medianArray, ciArray ] = ci_communityProperties( propsArray, alp )
%CI_COMMUNITYPROPERTIES Returns predictive intervals for community
%properties outlined in Hill et al. 2014
%
% Input
%   propsArray - cell array; each cell holds samples of a property;
%       if vector property (e.g. recurrence times, turnover times)
%           rows - state; cols - samples
%       if scalar property (e.g. damping ratio, convergence rate), 
%           cell is a row vector
%   alp - alpha level; default = 5%
% Output
%   muArray - means for each property
%   medianArray - medians 
%   ciArray - predictive intervals; default = 95% bounds

nprops = length(propsArray);
ciArray = {}; muArray = {}; medianArray = {}

if nargin < 2
    alp = 0.05;
end

for i = 1:nprops
    disp(['property',num2str(i)])
    tmp = propsArray{i};
    nr = size(tmp,1);
    tmp_ci = [];
    if nr > 1        
        for r = 1:nr
            tmp2 = tmp(r,:)';
            [a,b] = hpd_sim(tmp2, alp);
            tmp_ci(end+1,:) = [a,b];
        end
    else
        [a,b] = hpd_sim(tmp', alp);
        tmp_ci(1,1:2) = [a,b];
    end
    muArray{i} = mean(tmp,2);
    medianArray{i} = median(tmp,2);
    ciArray{i} = tmp_ci;
    
end

end

