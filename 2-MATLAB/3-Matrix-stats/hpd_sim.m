function [LB, UB] = hpd_sim(data, alpha)
%HPD_SIM Returns the bounds for highest posterior density region
% Purpose: 
% Compute the HPD region using random sample (posterior MCMC draws) 
% from some distribution
% -----------------------------------
% Algorithm: 
% Follows the method proposed by Chen and Shao (1998)
% It can only find the HPD of a single-modal distribution.
% -----------------------------------
% Arguments:
% data = random sample (posterior MCMC draws), a column vector 
% alpha = significance level, say 0.05, (In that case, 95% HPD interval will be computed)
% -----------------------------------
% Returns:
% LB = lower bound of the HPD interval
% UB = upper bound of the HPD interval
% 
% From the code written by Hang Qian, Iowa State University (matlabist@gmail.com)
% with minor modifications

if nargin < 2
    alpha = 0.05;
end

nobs = length(data);
cut = round(nobs * alpha);
span = nobs - cut;

data = sort(data);
[gap,ind] = min( data(span+1:nobs) - data(1:cut) );
LB = data(ind);
UB = data(ind + span);
