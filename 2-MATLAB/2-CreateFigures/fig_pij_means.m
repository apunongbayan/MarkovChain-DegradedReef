function [ output_args ] = fig_pijs_means( P_mu )
%FIG_PIJS_MEANS Summary of this function goes here
%   Detailed explanation goes here

figure();

%rescale data to be in the range 1-64
dmin = -7.08; dmax = 0;
rng = dmax - dmin;
d = 1 + 63*(log(P_mu)-dmin)/rng;    %scale the data
image(d);
caxis([2 64.0000]);
colormap(gray);
cmap = colormap;
cmap = flipud(cmap);                %flip color map
colormap(cmap);

ylabel({'fate'});
set(gca,...
    'YTickLabel',{'C: 1','MA: 2','ACA: 3','CAL: 4','OT: 5','HCJ: 6','HC: 7','SC: 8'},...
    'YTick',[1 2 3 4 5 6 7 8],...
    'YDir','reverse');

Ptmp = P_mu';
for ii = 1:8
    for jj = 1:8
        Pij = Ptmp(ii,jj);
        if Pij >= 0.10
            curcol = [1 1 1];
        else
            curcol = [0 0 0];
        end
        text('Position',[ii-0.35,jj],...
            'string',sprintf('%.3f',Pij),'fontsize',10,'Color',curcol);
    end
end
 
pos = get(gca,'position');
hc = colorbar('location','manual','position',...
     [pos(1)+pos(3)+.01 pos(2) .03 pos(4)]);
%where you want tick marks (unscaled data)
L = [0.001:.001:0.009, 0.01:.01:.09, 0.1:0.1:0.9, 1];
Lhat = {'0.001','','','','','','','','',...
    '0.01','','','','','','','','',...
    '0.1','','','','','','','','',...
    '1'};
l = 1 + 63*(log(L)-dmin)/rng;
set(hc,'Ytick',l,'YTicklabel',Lhat, 'ylim', [2 64.0000]);
xlabel({'state'});
title('hierarchical model');


end

