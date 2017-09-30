function [ output_args ] = fig_competitive_rank( R )
%FIG_COMPETITIVE_RANK Visualize the states that appear in rank 1, rank 2,
%etc.
%   Input: R - matrix; 
%       rows - states sorted in descending rank order; 
%       cols - samples

[ncompetitors,nsamps] = size(R);

xvals = unique(R);
xbins = min(xvals):1:max(xvals);

figure()
labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)'}
for i = 1:ncompetitors
    subplot(ncompetitors,1,i)
    [freqs,xcenters] = hist(R(i,:),xbins)
    px = freqs / sum(freqs);
    bar(xbins,px,'FaceColor',[0.5 0.5 0.5],...
            'EdgeColor',[0.5 0.5 0.5])
    xlim([1,9]);
    text = [labels{i},' Rank ',num2str(i)];
    title(text)
    if i==4
        ylabel('p(x)','Fontsize',12)
    end
end
xlabel('state','Fontsize',12)

end

