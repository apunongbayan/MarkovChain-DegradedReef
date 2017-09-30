function [ output_args ] = addtrajectories( treatnum, statenum)
%ADDTRAJECTORIES Adds observed trajectories in treatment plots to an
%existing figure
%   

load matrices-2015-01-05.mat 
cover = {obsCover_n_all, obsCover_e, obsCover_h, obsCover_hs};

treats = {'n','e','e+hc','e+hc+sc'}
%ym = [1 .8 .4 1 .1 .02 .3 .2]
ym = [1 .8 .5 1 .1 .02 .4 .2]
%ym = [1 1 1 1 1 1 1 1]
sts = {'E','MA','ACA','CAL','OT','HCJ','HC','SC'};
marks = {'s','o','^','v'}
msizes = [7 7 7 7]


ci = cover{treatnum};
nq = length(ci)
cc=bone(nq+3)

if nargin < 2
    statenum = 7;
end

% individual plots
% for i = 1:nq;
%     ctmp = ci{i}
%     if length(statenum) > 1
%         plot(ctmp(2:end,4)-1, sum(ctmp(2:end,statenum+4),2)*100,...
%             'Marker',marks{treatnum},'MarkerEdgeColor',cc(i+2,:),...
%             'MarkerFaceColor','none','MarkerSize',msizes(treatnum),...
%             'LineStyle','none','LineWidth',1,'Color',cc(i,:)); hold on; 
%     else        
%         plot(ctmp(2:end,4)-1, ctmp(2:end,4+statenum)*100,...
%             'Marker',marks{treatnum},'MarkerEdgeColor',cc(i+1,:),...
%             'MarkerFaceColor','none','MarkerSize',msizes(treatnum),...
%             'LineStyle',':','LineWidth',1,'Color',cc(i,:)); hold on; 
%     end
% end


% means and se
xi = obsCover_treatmu{treatnum}; sei = obsCover_treatse{treatnum};
tvec = xi(2:end,3)-1
% errorbar(tvec,xi(2:end,3+statenum)*100, sei(2:end,3+statenum)*200,...
%             'Marker',marks{treatnum},'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[1 1 1],'MarkerSize',msizes(treatnum),...
%             'LineStyle','none','Color',cc(treatnum,:));
plot(tvec,sum(xi(2:end,3+statenum),2)*100,...
            'Marker',marks{treatnum},'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],'MarkerSize',msizes(treatnum),...
            'LineStyle','none','Color',cc(treatnum,:));


% model equilibrium cover - be sure to change when using another model 
% statvec = statdist(jags_p{2});
% plot(tvec,repmat(statvec(statenum),1,length(tvec))*100,'k--')

xlim([0,44])


end

