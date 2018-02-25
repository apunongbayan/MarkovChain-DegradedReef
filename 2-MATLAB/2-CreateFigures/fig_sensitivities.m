function [ ] = fig_sensitivities( SensitivityMatrices )
%FIG_SENSITIVITIES Create figure of sensitivities and their envelopes
%   Input
%       SensitivityMatrices
%           (x,y) - sensitivity values
%           z - sensitivity matrix specific to a transition probability
%           matrix
%

% Consider only the real parts (since imaginary parts are negligible (=0.0000i)
Sreal = real(SensitivityMatrices);
% Get median values for sensitivity to each transition probability pij
Smedian = median(Sreal, 3);

[nx,ny,nz] = size(Sreal);

% Calculate highest posterior density intervals
UB95 = zeros(nx,nx); LB95 = zeros(nx,nx);
UB50 = zeros(nx,nx); LB50 = zeros(nx,nx);
alph95 = 0.05;    % 95% HPD
alph50 = 0.50;    % 50% HPD

for j = 1:8
    for i = 1:8
        tmp(:,1) = Sreal(i,j,:);    % sampled values as column vector
        [LB95(i,j), UB95(i,j)] = hpd_sim(tmp, alph95); % calculate 95% HPD                
        [LB50(i,j), UB50(i,j)] = hpd_sim(tmp, alph50); % calculate 50% HPD
    end
end

%Convert to vectors
svec = vec(Smedian);    % median values
uvec95 = vec(UB95);     % 95% HPD
lvec95 = vec(LB95);
uvec50 = vec(UB50);     % 50% HPD
lvec50 = vec(LB50);

labelvec = {'11','21','31','41','51','61','71','81',...
    '12','22','32','42','52','62','72','82',...
    '13','23','33','43','53','63','73','83',...
    '14','24','34','44','54','64','74','84',...
    '15','25','35','45','55','65','75','85',...
    '16','26','36','46','56','66','76','86',...
    '17','27','37','47','57','67','77','87',...
    '18','28','38','48','58','68','78','88',};

% sort mean sensitivities according to decreasing absolute magnitude, 
% tmp = abs(svec);
% [sens,idum] = sort(tmp,'descend');
idum = 1:length(svec);

% find sensitivities that are consistent in sign
% (i.e. 95% HPD do not overlap both + and - values)
consistency_check = lvec95(idum).*uvec95(idum);
%consistent = idum(find(consistency_check>0));
thresh = 0.02;
consistent = idum(find(consistency_check>0 & abs(svec) > thresh))
inconsistent = idum(find(consistency_check<=0 | abs(svec) <= thresh));


% Plot consistent sensitivities
figure()                        % used for manuscript figure
topid = consistent;     
xi = length(topid);
plot([0,xi+1],[0,0], 'k--'); hold on;   % plot line y = 0 
                    
for i = 1:xi
    plot([i i],[lvec95(topid(i)) uvec95(topid(i))],'k-','LineWidth',1)
    hold on;
    plot([i i],[lvec50(topid(i)) uvec50(topid(i))],'k-','LineWidth',2)
end
plot(1:xi,svec(topid),'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1]);

% % check
% labelvec(topid)
% svec(topid)
% lvec95(topid)
% uvec95(topid)
% min(abs(svec(topid)))

% labels = cellstr(num2str(topid));  %' # labels correspond to their order
oddxi = 1:2:length(topid);
evenxi = 2:2:length(topid);
text(oddxi, uvec95(topid(oddxi))+0.04, labelvec(topid(oddxi)), 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',10)
text(evenxi, uvec95(topid(evenxi))+0.07, labelvec(topid(evenxi)), 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',10)

xlabel('Transition (p_i_j)','FontSize',12)
ylabel('Sensitivity of stationary coral cover (s_i_j)','FontSize',12)
set(gca,'XTickLabel',{'',''},'XTick',[0 length(topid)+1])
xlim([0 length(topid)+1])
ylim([-0.7 2.8])


% Supplemental figure: small and/or inconsistent sensitivities 
figure()
% topid = [consistent(24:end,1); inconsistent(:,1)] 
topid = inconsistent;
xi = length(topid);
plot([0,xi+1],[0,0], 'k--'); hold on; 

for i = 1:xi
    plot([i i],[lvec95(topid(i)) uvec95(topid(i))],'k-','LineWidth',1)
    hold on;
    plot([i i],[lvec50(topid(i)) uvec50(topid(i))],'k-','LineWidth',2)
end
plot(1:xi,svec(topid),'ko','MarkerSize',6,'MarkerFaceColor',[1 1 1]);

% labels = cellstr(num2str(topid));  %' # labels correspond to their order
% labels = labelvec(topid)
% text(1:xi, svec(topid), labels, 'VerticalAlignment','bottom', ...
%     'HorizontalAlignment','right','FontSize',12)
oddxi = 1:2:length(topid);
evenxi = 2:2:length(topid);
text(oddxi, uvec95(topid(oddxi))+0.01, labelvec(topid(oddxi)), 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',10)
text(evenxi, uvec95(topid(evenxi))+0.02, labelvec(topid(evenxi)), 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',10)

xlabel('Transition (p_i_j)','FontSize',12)
ylabel('Sensitivity of stationary coral cover (s_i_j)','FontSize',12)
set(gca,'XTickLabel',{'',''},'XTick',[0 length(topid)+1])
xlim([0 length(topid)+1])
ylim([-0.3 0.3])

end

