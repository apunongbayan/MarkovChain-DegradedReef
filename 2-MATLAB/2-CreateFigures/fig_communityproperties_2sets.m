function [ output_args ] = fig_communityproperties_2sets( samplesArray1,samplesArray2 )
%FIG_COMMUNITYPROPERTIES Creates a figure for a subset of properties
%   Means and prediction intervals (50%, 95%)
%   Input
%       sampleArray - each cell contains samples of a given property/index

% Set 1
% 50% HPD intervals
[mu_cprops1, median_cprops1, ci50_cprops1] = ci_communityProperties(samplesArray1,0.5);
% 95% HPD intervals
[mu_cprops1, median_cprops1, ci_cprops1] = ci_communityProperties(samplesArray1,0.05);

% Set 2
% 50% HPD intervals
[mu_cprops2, median_cprops2, ci50_cprops2] = ci_communityProperties(samplesArray2,0.5);
% 95% HPD intervals
[mu_cprops2, median_cprops2, ci_cprops2] = ci_communityProperties(samplesArray2,0.05);




figure();
subplot1(2, 2, 'Gap', [.08 .03],'YTickL','All')
dx = 0.18;

% % Turnover times
pnum = 3; 
subplot1(1)
% Set 1
for i=1:8    
    plot([i-dx,i-dx], ci_cprops1{pnum}(i,:)*12, 'k-','LineWidth',1.5); hold on
    %plot([i-dx,i-dx], ci50_cprops1{pnum}(i,:)*12, 'k-','LineWidth',3); hold on
    plot(i-dx,mu_cprops1{pnum}(i)*12,'ko','MarkerSize',6,'MarkerFaceColor',[0 0 0])
    %plot(i-dx,median_cprops1{pnum}(i)*12,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
% Set 2
for i=1:8
    plot([i+dx,i+dx], ci_cprops2{pnum}(i,:)*12, '-','LineWidth',1.5,...
        'Color',[0.65 0.65 0.65]); hold on
    %plot([i+dx,i+dx], ci50_cprops2{pnum}(i,:)*12, 'k-','LineWidth',3); hold on
    plot(i+dx,mu_cprops2{pnum}(i)*12,'o','MarkerSize',6,...
        'MarkerFaceColor',[0.65 0.65 0.65],...
        'MarkerEdgeColor',[0.65 0.65 0.65]);
    %plot(i+dx,median_cprops2{pnum}(i)*12,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
xlim([0,9])
ylabel('Turnover time (months)','FontSize',12);
set(gca,'XTickLabel',{'','','','','','','','','',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])
text(.10,.93,'(a)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);


% % Recurrence times
pnum = 2; 
subplot1(2)
% Set 1
for i=1:8    
    plot([i-dx,i-dx], ci_cprops1{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i-dx,i-dx], ci50_cprops1{pnum}(i,:), 'k-','LineWidth',3); 
    plot(i-dx,mu_cprops1{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[0 0 0]);
    %plot(i-dx,median_cprops1{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
% Set 2
for i=1:8
    plot([i+dx,i+dx], ci_cprops2{pnum}(i,:), '-','LineWidth',1.5,...
        'Color',[0.65 0.65 0.65]); hold on
    %plot([i+dx,i+dx], ci50_cprops2{pnum}(i,:), 'k-','LineWidth',3); 
    plot(i+dx,mu_cprops2{pnum}(i),'o','MarkerSize',6,...
        'MarkerFaceColor',[0.65 0.65 0.65],...
        'MarkerEdgeColor',[0.65 0.65 0.65]);
    %plot(i+dx,median_cprops2{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
xlim([0,9])
ylabel('Recurrence time (years)','FontSize',12);
text(.10,.93,'(b)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
set(gca,'YScale','log','YMinorTick','on',...
    'XTickLabel',{'','','','','','','','','',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])


% % Entropy of transitions given current state
pnum = 12;
subplot1(3)
%average entropy over states
plot([0:9], repmat(mu_cprops1{5},[1 10]), 'k:'); hold on     
plot([0:9], repmat(mu_cprops2{5},[1 10]), ':','Color',[0.65 0.65 0.65]); 
% set 1
for i=1:8    
    plot([i-dx,i-dx], ci_cprops1{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i-dx,i-dx], ci50_cprops1{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i-dx,mu_cprops1{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[0 0 0]);
    %plot(i-dx,median_cprops1{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
% set 2
for i=1:8
    plot([i+dx,i+dx], ci_cprops2{pnum}(i,:),'-','LineWidth',1.5,...
        'Color',[0.65 0.65 0.65]); hold on
    %plot([i+dx,i+dx], ci50_cprops2{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i+dx,mu_cprops2{pnum}(i),'o','MarkerSize',6,...
        'MarkerFaceColor',[0.65 0.65 0.65],...
        'MarkerEdgeColor',[0.65 0.65 0.65]);
    %plot(i+dx,median_cprops2{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
text(.10,.93,'(c)','Units', 'Normalized',...
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
% set 1
for i=2:8
    plot([i-dx,i-dx], ci_cprops1{pnum}(i,:), 'k-','LineWidth',1.5); hold on
    %plot([i-dx,i-dx], ci50_cprops1{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i-dx,mu_cprops1{pnum}(i),'ko','MarkerSize',6,'MarkerFaceColor',[0 0 0])
    %plot(i-dx,median_cprops1{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
% set 2
for i=2:8
    plot([i+dx,i+dx], ci_cprops2{pnum}(i,:), '-','LineWidth',1.5,...
        'Color',[0.65 0.65 0.65]); hold on
    %plot([i+dx,i+dx], ci50_cprops2{pnum}(i,:), 'k-','LineWidth',3); hold on
    plot(i+dx,mu_cprops2{pnum}(i),'o','MarkerSize',6,...
        'MarkerFaceColor',[0.65 0.65 0.65],...
        'MarkerEdgeColor',[0.65 0.65 0.65]);
    %plot(i+dx,median_cprops2{pnum}(i),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1]);
end
text(.10,.93,'(d)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
ylim([0,0.2])
xlim([0,9])
xlabel('State','FontSize',12)
ylabel('Replacement by','FontSize',12);
set(gca,'XTickLabel',{'','C','MA','ACA','CAL','OT','HCJ','HC','SC',''},...
    'XTick',[0 1 2 3 4 5 6 7 8 9])


end

