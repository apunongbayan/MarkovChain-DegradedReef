function [  ] = fig_treatmentpreds( trajmx, tvec )
%FIG_TREATMENTPREDS Plot prediction envelope for a given treatment
%   Subplots pertaining to each state
% Input
%   trajmx - 3d matrix, x - replicate orbit, y - time, z - state;

[ns, nt, nr] = size(trajmx);

figure(); hold on;
subplot1(4, 2, 'Gap', [.02 .02],'YTickL','Margin','XTickL','Margin');
%tvec = 0:1:nt-1;
flabs = {'(a) C','(b) MA','(c) ACA','(d) CAL',...
    '(e) OT','(f) HCJ','(g) HC','(h) SC'};
yub = [100 100 100 100 20 20 20 20];

ppimed = [];
ppimed(:,1) = tvec;

trajs = zeros(nr,nt,ns);
lbs = zeros(nt,ns);

for si = 1:ns
    for ri = 1:nr
        trajs(ri,:,si) = trajmx(si,:,ri);    
    end
    
    % calc 95% HPD
    statecover = trajs(:,:,si);
    alpha = .05;
    lb = [];
    ub = [];
    for ti = 1:nt
        [lb(end+1), ub(end+1)] = hpd_sim(statecover(:,ti),alpha);          
    end
    lbs(:,si) = lb; 
    ubs(:,si) = ub; 
    subplot1(si); hold on;
    
    % Plot 95% HPD for trajectory
    plot(tvec,lb*100,'k--'); hold on; plot(tvec,ub*100,'k--');
    
    %Option 1 - median
    %ppimed(:,end+1) = median(statecover)*100;
    %plot(tvec,median(statecover)*100,'k-');    
    
    %Option 2 - mean
    ppimed(:,end+1) = mean(statecover)*100;
    plot(tvec,mean(statecover)*100,'k-'); 
    
    % Add observed trajectory in treatment 2 plots
%     addtrajectories(2,si);
    
    % Figure related stuff
    text(.65,.95,[flabs{si}],'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    ylim([0 yub(si)]);
    if ismember(si,[1 3 5 7])
        ylabel('Cover (%)', 'FontSize', 11)           
    end
    
    if ismember(si,[7 8])
        xlabel('Time (months)', 'FontSize', 11)           
    end

end

% filename = ['unkAlpha_mean_treat2.mat'];
% save(filename,'ppimed','lbs','ubs');

end

