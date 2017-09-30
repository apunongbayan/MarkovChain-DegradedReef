function [ output_args ] = fig_paramposterior_bars( samples, polycol,...
    xlimon, pind, lstyle)
%FIG_PARAMPOSTERIOR_LINES Create histograms for densities of probabilities; 
%already includes label for each subplot
%for concentration parameters use fig_paramposterior_lines2

% Input
% samples - matrix, row - sample , col - param
% polycol - color of bar
% xlimon - indicator for setting x limits to maximum of observed
% pind - indicator for samples that are constrained to [0,1]
% lstyle - line style


samples = rearrange_jagssamples(samples);

ts = {'', '', '', '', '', '', '', '',...
    '', 'MA', 'ACA', 'CAL', 'OT', 'HCJ', 'HC', 'SC',...
    '', '', '', '', '', '', '', '',...
    '', '', '', '', '', '', '', '',...
    '', '', '', '', '', '', '', '',...
    '', '', '', '', '', '', '', '',...
    '', '', '', '', '', '', '', '',...
    '', '', '', '', '', '', '', ''};

flabs = {'_1','_2','_3','_4','_5','_6','_7','_8'};

if nargin < 5
    lstyle = '-';
end

[xx,yy] = size(samples);
c = ceil(sqrt(yy));

if pind == 1            %WHEN WORKING WITH PROBABILITIES USE THIS
    subplot1(c, c, 'Gap', [.035 .02],'YTickL','All')
    sometext = '{\itI}';
else                    %WHEN WORKING WITH ALPHAS USE THIS
    subplot1(c, c, 'Gap', [0.035 .02], 'XTickL', 'All', 'YTickL','All')
    sometext = 'a';
end


for yi=1:yy
    %subplot(c,c,yi);
    subplot1(yi);
    theta_x = mean(samples(:,yi));
    bwid = theta_x * 0.1;
    xbins = .00:.1:1;
    if pind == 1
        [freqs, xcenters ] = hist(samples(:,yi),xbins);
    else        
        [freqs, xcenters ] = hist(samples(:,yi),0:bwid:2000);  
    end
    
    dx = xcenters(2) - xcenters(1);
    fi = freqs * dx;
    ntot = sum(fi);
    px = (fi / ntot) / dx;   %probability density 

    bar(xcenters,px,'FaceColor',polycol,...
            'EdgeColor',polycol); 
    hold on;
    plot(median(samples(:,yi)),0,'ko',...
        'MarkerSize',7,'MarkerFaceColor','r')
        
    %Used linewidth of 2 for the paper figures; 1.5 for supplementary
    %figures
        
    if xlimon == 1
        idum = find(freqs > 0);
        xlim([0 xcenters(idum(end))+dx/2]);
        %ylim([-0.01 max(px)]);  
        ylim([-0.00001 max(px)]);
    end

    si = mod(yi,8) ;
        if si == 0
            si = 8;
        end            
    fi = floor((yi-1)/8)+1; 
    text(.55,.95,[sometext,flabs{fi},flabs{si}],'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',14, 'FontName','Times New Roman');

    if yi == 33
        ylabel(' ','FontSize',14);
    end   


    %WHEN WORKING WITH PROBABILITIES
    if pind==1
        xlim([-0.025 1])
        ylim([-0.01 max(px)]);  
    end
    
    if and(ismember(yi,[57:64]),pind==1)
        xlabel('')
        set(gca,'XTickLabel',{'0','0.5','1'},'XTick',[0 0.5 1])
    end
    
    
    if ismember(yi,[9:16])
        title(ts{yi},'FontSize',14)
    end
    
    if yi == 60
        xlabel(' ','FontSize',14);
    end 
    hold on;
    
end

% h = text(1,7,'Posterior density');
% set(h, 'Rotation', 90);

% 13cm width? 
end
