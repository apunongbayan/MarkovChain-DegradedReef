function [ rsamples ] = rearrange_jagssamples( samples )
%REARRANGE_JAGSSAMPLES Takes jags samples and rearranges columns
% Rearranges columns so that histograms could be plotted with 
% columns containing the state, and rows the fate
% Used by: fig_paramposterior_lines.m

nk = 8;  %8 state model
nparms = 64;

[nr,nc] = size(samples);
rsamples = [];

for i=[1:nk, 0]
    for j = 1:nparms
        if mod(j,8) == i
            rsamples(:,end+1) = samples(:,j);
        end
    end
end
            
            


end

