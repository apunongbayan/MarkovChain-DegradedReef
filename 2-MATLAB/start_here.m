%% The following script is for the simulations and matrix analyses

% Fixed effects model with uniform Dirichlet prior - P samples

close all; clear all;

% Import MCMC samples for the parameters (the output from R)
% columns 1-64 are transition probabilities pij's
% columns 65-72 are initial state pi_0
% Make sure that Input folder  has been added to path
filename = 'jags-output-t0-26-v1-fixed.csv';
header_row = 1;
posterior_draws = csvread(filename,header_row);

PmatSamplesRows = posterior_draws(:,1:64);      %transition probabilities
InitStateDistRows = posterior_draws(:,65:end);  %initial state distribution

file = 'FixedEffect_P_samples.mat';
save(file,'PmatSamplesRows','InitStateDistRows');

%% Categorical-Dirichlet model for patch state transitions
% alpha samples

close all; clear all;

% Import MCMC samples for the parameters (the output from R)
% columns 1-64 are Dirichlet priors for transition probabilities pij's
% columns 65-72 are Dirichlet priors for initial state pi_0
% Make sure that Input folder  has been added to path
filename = 'jags-output-t0-26-v2-multidirch.csv';
header_row = 1;
A = csvread(filename,header_row);

% MCMC samples (cols - params; rows - samples)
Ap_samples = A(:,1:64);     % Dirichlet priors for pijs
Ax_samples = A(:,65:end);   % Dirichlet priors for initial state distribution 

% Calculate the posterior means
A_mu = mean(A,1);           
Ap_mu = vec2mat(A_mu(1:64),8)'; 
Ax_mu = A_mu(65:end);

% Expected values
P_mu = Ap_mu ./ repmat(sum(Ap_mu),8,1);     % transition probability matrix
x_mu = Ax_mu / sum(Ax_mu);                  % initial state               

% Save variables for later use
file2 = ['CatDirch_params.mat'];
save(file2,'Ap_samples','Ax_samples','Ap_mu','Ax_mu','P_mu','x_mu');

%% Inspect mean pijs (Fig. 2 subplot)

close all; clear all;
load CatDirch_params.mat 

% Inspect relative magnitudes of the mean transition probabilities
fig_pij_means(P_mu);


%% Sample transition probabilities
% Option 1 - given uncertainty about Dirichlet parameters

close all; clear all; 
load CatDirch_params.mat 

[nr,nc] = size(Ap_samples);
ns = sqrt(nc);

nsets = 30000; % Number of Dirichlet parameter sets to sample
nq = 3;     % number of replicate quadrats given Dirichlet parameter set
PmatSamplesRows = [];
InitStateDistRows = [];

for n = 1:nsets
    ri = randi(nr,1);
    A = vec2mat(Ap_samples(ri,:),ns)';
    a0 = Ax_samples(ri,:);
    for q = 1:nq
        % sample transition prob. matrix
        P_i = samplePmatfromAlpha(A);
        % convert to row and save to list
        PmatSamplesRows(end+1,:) = P_i(:)';
        
        %sample initial state distribution
        InitStateDistRows(end+1,:) = sampleInitDistfromAlpha(a0);
    end
end

file = 'UnkAlpha_P_samples.mat';
save(file,'PmatSamplesRows','InitStateDistRows');

%% Sample transition probabilities 
% Option 2 - using the posterior means of Dirichlet parameters

close all; clear all;
load CatDirch_params.mat 

nq = 90000;    % number of replicate quadrats given Dirichlet parameter set
PmatSamplesRows = [];
InitStateDistRows = [];

for n = 1:nq
    P_i = samplePmatfromAlpha(Ap_mu);
    PmatSamplesRows(end+1,:) = P_i(:)';   
    InitStateDistRows(end+1,:) = sampleInitDistfromAlpha(Ax_mu);
end

file = 'MeanAlpha_P_samples.mat';
save(file,'PmatSamplesRows','InitStateDistRows');

%% Inspect the posterior densities of pijs 

close all; clear all; 

% % Option 0 - fixed effects model
% load FixedEffect_P_samples.mat

% Option 1 - given parameter uncertainty 
load UnkAlpha_P_samples.mat

% % Option 2 - when drawn from posterior means of Dirichlet parameters
% load MeanAlpha_P_samples.mat

linecolor = [0 0 0];
pij_indicator = 1;
xaxis_maxval = 0;

fig_paramposterior_lines(PmatSamplesRows,linecolor, xaxis_maxval, ...
    pij_indicator, '-');

%% Inspect the posterior densities of Dirichlet parameters

close all; clear all;
load CatDirch_params.mat 

linecolor = [0 0 0];
pij_indicator = 0;
xaxis_maxval = 1;

fig_paramposterior_lines(Ap_samples,linecolor, xaxis_maxval, ...
    pij_indicator, '-');

%% Calculate the stationary distribution - Fig. 3 

close all; clear all; 

% % Option 0: Fixed effects model
% load FixedEffect_P_samples.mat
% file = 'FixedEffect_w1_samples.mat';  % output file

% Option 1: Given parameter uncertainty
load UnkAlpha_P_samples.mat 
file = 'UnkAlpha_w1_samples.mat';   %output file

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_P_samples.mat
% file = 'MeanAlpha_w1_samples.mat'; %output file

[nx,ny] = size(PmatSamplesRows);
ns = sqrt(ny);

w1_samples = zeros(ns,nx);

for i = 1:nx
    w1_samples(:,i) = statdist(reshape(PmatSamplesRows(i,:),ns,ns));
end

fig_histstatdist(w1_samples');
save(file,'w1_samples');

%% Projected trajectories in natural / unmanipulated plots

clear all; close all;

% % Option 0:  Fixed effects model
% load FixedEffect_P_samples.mat 
% file = 'FixedEffect_projections_naturalplots.mat';

% % Option 1: Given uncertainty about Dirichlet parameters
load UnkAlpha_P_samples.mat 
file = 'UnkAlpha_projections_naturalplots.mat';

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_P_samples.mat
% file = 'MeanAlpha_projections_naturalplots.mat';

time_end = 26;                  % months

[trajs, tvec] = stateTrajectoriesRandInit(PmatSamplesRows,...
    InitStateDistRows,time_end);

save(file,'trajs','tvec');

fig_treatmentpreds(trajs, tvec)

%% Projected trajectories in empty plots

clear all; close all;

% % Option 0:  Fixed effects model
% load FixedEffect_P_samples.mat 
% file = 'FixedEffect_projections_clearedplots.mat';

% % Option 1: Given uncertainty about Dirichlet parameters
load UnkAlpha_P_samples.mat 
file = 'UnkAlpha_projections_clearedplots.mat';

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_P_samples.mat
% file = 'MeanAlpha_projections_clearedplots.mat';

init_x = [1 0 0 0 0 0 0 0]';    % initial state distribution; 100% C
time_end = 44;                  % months

[trajs, tvec] = stateTrajectories(PmatSamplesRows,init_x,time_end);

save(file,'trajs','tvec');

%% Create figure of predicted trajectories starting from 100% crustose algae

clear all; close all;

% % Option 0:  Fixed effects model
% load FixedEffect_projections_clearedplots.mat

% % Option 1 - unknown Dirichlet parameters
load UnkAlpha_projections_clearedplots.mat

% % Option 2 - posterior mean Dirichlet parameters
% load MeanAlpha_projections_clearedplots.mat

fig_treatmentpreds(trajs, tvec)

%% Calculate successional indices (means and 95% predictive intervals)

close all; clear all; 

% % Option 0: Fixed effects
% load FixedEffect_P_samples.mat
% file = 'FixedEffect_indices.mat';

% Option 1: Given uncertainty about Dirichlet parameters
load UnkAlpha_P_samples.mat 
file = 'UnkAlpha_indices.mat';    % output file

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_P_samples.mat
% file = 'MeanAlpha_indices.mat';     % output file     

interval_length = 1;    % resolution of P (in months)
empty_state_index = 1;  % CCA or empty state

% Calculate successional indices (recurrence times, turnover times, etc.
% for all samples of transition probability matrices
[all_indices_samples,propertyNames] = ...
    communityProperties(PmatSamplesRows,interval_length,empty_state_index);
% Calculate the highest posterior density regions (95% and 50% HPDs)
[means,medians,hpd95] = ci_communityProperties(all_indices_samples);
[means,medians,hpd50] = ci_communityProperties(all_indices_samples,0.5);

% Save samples and stats to file
save(file,'propertyNames','all_indices_samples','means','medians',...
    'hpd95','hpd50');


%% Plot a subset of successional indices
% Just 1 model

close all; clear all; 

% % Option 1: Given uncertainty about Dirichlet parameters
load UnkAlpha_indices.mat

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_indices.mat

fig_communityproperties(all_indices_samples)


%% Plot a subset of successional indices
% Latest version: 2 sets (with or without parameter uncertainty)

% % close all; clear all; 
% % 
% % % % Option 1: Given uncertainty about Dirichlet parameters
% % load UnkAlpha_indices.mat
% % props1 = all_indices_samples;
% % 
% % % Option 2: Given posterior means for Dirichlet parameters
% % load MeanAlpha_indices.mat
% % props2 = all_indices_samples;
% % 
% % fig_communityproperties_2sets(props1,props2)


%% Sensitivity analysis

close all; clear all;

% Option 1: Given uncertainty about Dirichlet parameters
load UnkAlpha_P_samples.mat 
file = 'UnkAlpha_sensmatsamples.mat';

% % Option 2: Given posterior means for Dirichlet parameters
% load MeanAlpha_P_samples.mat
% file = 'MeanAlpha_sensmatsamples.mat';

% Calculate sensitivity uncertainty using smaller number of samples
PmatSamplesSub = PmatSamplesRows(1:10000,:);

% pick only the sensitivities of hard coral states (juveniles and adults)
weight = [0 0 0 0 0 1 1 0]; 
% calculate sensitivity matrix for each sampled transition probability mx.
Smats = sensitivities(PmatSamplesSub,weight);

save(file,'Smats');

%% Create figure for sensitivities

close all; clear all;

% Option 1 - accounts for parameter uncertainty
load UnkAlpha_sensmatsamples.mat

% % Option 2 - using posterior means of parameters
% load MeanAlpha_sensmatsamples.mat

fig_sensitivities(Smats)




