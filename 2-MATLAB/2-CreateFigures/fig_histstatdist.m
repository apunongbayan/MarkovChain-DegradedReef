function [ output_args ] = fig_histstatdist( samples )
%FIG_HISTSTATDIST Creates histogram plots for parameters in mx. samples
%
% Input
%  samples - a matrix; col - parameter; row - samples


xmaxs = [0.7, 0.3, 0.4, 0.7, 0.15, 0.1, 0.15, 0.1]*100;
flabs = {'(a) C','(b) MA','(c) ACA','(d) CAL','(e) OT',...
    '(f) HCJ','(g) HC','(h) SC'};

figure('Renderer','painters','Color',[1 1 1]);
subplot1(2, 4, 'Gap', [.04 .06], 'XTickL','All','YTickL','All'); %gap 0.4,0.6
for i = 1:8
    subplot1(i);

    %stationary dist histogram
    %hist(samples(:,i)*100,100);
    [freqs,xcenters] = hist(samples(:,i)*100,100);
    bar(xcenters,freqs,'FaceColor',[0.6 0.6 0.6],...
            'EdgeColor',[0.6 0.6 0.6]); hold on
       
    hold on;
    %h = findobj(gca,'Type','patch');
    %set(h,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]);
   
    % central measures
    mm = mean(samples(:,i)*100);    
    [di, mmx] = min(abs(xcenters - (mm/100)));
    
    % shortest interval; HPD
    mo = median(samples(:,i)*100);
    [di, mox] = min(abs(xcenters - (mo/100)));
    %mo = mode(samples(:,i));
    [lb, ub] = hpd_sim(samples(:,i)*100,0.05);
    [di, lbx] = min(abs(xcenters - (lb/100)));
    [di, ubx] = min(abs(xcenters - (ub/100)));

%     %old version - 95%CI as vertical lines
    plot([lb, lb],[0,max(freqs)],'k--', 'LineWidth',1);
    plot([ub, ub],[0,max(freqs)],'k--', 'LineWidth',1);
    plot([mo, mo],[0,max(freqs)],'Color',[0.1 0.1 0.1],'LineStyle','-',...
        'LineWidth',2);
    
    %label
    if mod(i,4) == 1
        ylabel('Density','FontSize',11);
    end
    if i == 5 | i == 6 |i == 7 | i == 8
        xlabel('Cover (%)','FontSize',11);
    end
    
    xlim([0 xmaxs(i)]);
    
    % before .7,.95
    text(.65,.97,flabs{i},'Units', 'Normalized',...
        'VerticalAlignment', 'Top','FontSize',11);
    

    
    % force scientific notation on axes, for uniformity
    % code from: https://www.mathworks.com/matlabcentral/answers/158707-force-scientific-notation-in-axes
    % get ylim
    yl=ylim;
    
    % get order of magnitude
    e=log10(yl(2));
    e=sign(e)*floor(abs(e));
    
    % get and rescale yticks
    yt=get(gca,'ytick')/10^e;
    
    % create tick labels
    ytl=cell(size(yt));
    for j=1:length(yt)
        % the space after the percent gives the same size to positive and
        % negative numbers. The number of decimal digits can be changed.
        ytl{j}=sprintf('% 1.1f',yt(j));
    end
    
    % set tick labels
    set(gca,'yticklabel',ytl);
    
    % place order of magnitude
    fs = get(gca,'fontsize');
    set(gca,'units','normalized');
    xl = xlim;
    %text(xl(1),yl(2),sprintf('x 10^%d',e),...
    %text(xl(1),yl(2),sprintf('\\times10^%d',e),...
    text(xl(1),yl(2),sprintf('\\times10^{%d}',e),...
        'fontsize',fs,'VerticalAlignment','bottom');
    
end

%tightfig()


end

