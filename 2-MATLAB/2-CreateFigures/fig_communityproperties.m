function [ output_args ] = fig_communityproperties( samplesArray )
%FIG_COMMUNITYPROPERTIES Creates a figure for a subset of properties
%   Means and prediction intervals (50%, 95%)
%   Input
%       sampleArray - each cell contains samples of a given property/index


% 50% HPD intervals
[mu_cprops, median_cprops, ci50_cprops] = ci_communityProperties(samplesArray,0.5);
% 95% HPD intervals
[mu_cprops, median_cprops, ci_cprops] = ci_communityProperties(samplesArray,0.05);

figure();
subplot1(2, 2, 'Gap', [.08 .03],'YTickL','All')

% % Turnover times
pnum = 3; 
subplot1(1)
for i=1:8    
    plot([i,i], ci_cprops{pnum}(i,:)*12, 'k-','LineWidth',1.5); hold on
    %plot([i,i], ci50_cprops{pnum}(i,:)*12, 'k-','LineWidth',3); hold on
    plot(i,mu_cprops{pnum}(i)*12,'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1])
    %plot(i,median_cprops{pnum}(i)*12,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
xlim([0,9])
ylabel('Turnover time (months)','FontSize',12);
set(gca,'XTickLabel',{'','','','','','','','','',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])
text(.10,.95,'(a)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);


% % Recurrence times
pnum = 2; 
subplot1(2)
for i=1:8    
    plot([i,i], ci_cprops{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i,i], ci50_cprops{pnum}(i,:), 'k-','LineWidth',3); 
    plot(i,mu_cprops{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1]);
    %plot(i,median_cprops{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
xlim([0,9])
ylabel('Recurrence time (years)','FontSize',12);
text(.10,.95,'(b)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
set(gca,'YScale','log','YMinorTick','on',...
    'XTickLabel',{'','','','','','','','','',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])


% % Entropy of transitions given current state
pnum = 12;
subplot1(3)
for i=1:8    
    plot([i,i], ci_cprops{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i,i], ci50_cprops{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i,mu_cprops{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1]);
    %plot(i,median_cprops{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
plot([0:9], repmat(mu_cprops{5},[1 10]), 'k--')     %average entropy over states
text(.10,.95,'(c)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
xlim([0,9])
ylim([0.2,0.8])
xlabel('State','FontSize',12)
ylabel('Entropy','FontSize',12)
set(gca,'XTickLabel',{'','C','MA','ACA','CAL','OT','HCJ','HC','SC',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])

        
% % P(Replacement by j) averaged over all competing states (i != j)
pnum = 11;
subplot1(4)
for i=2:8    
    plot([i,i], ci_cprops{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i,i], ci50_cprops{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i,mu_cprops{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1])
    %plot(i,median_cprops{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
text(.10,.95,'(d)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
ylim([0,0.2])
xlim([0,9])
xlabel('State','FontSize',12)
ylabel('Replacement by','FontSize',12);
set(gca,'XTickLabel',{'','C','MA','ACA','CAL','OT','HCJ','HC','SC',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])




end

