% DDM for ONE subject ONE parameters
%%

    % Install wiener distribution package
    % Define priors (from paper)
    % Generate posterior for a few individuals look at the comparison of
    % distributions 

    % Recode reaction time as + or - for pump or cashout --> x is
    % choice-coded RT
    % Weiner diffusion model using jagssmples
    
    % x[t,p] ~ dwiener(alpha[t,p], tau[t,p], beta[t,p], delta[t,p])
    % Where t is trial index and p is pump index
    
    % To generatve new samples
    % px[t,p] ~ dwiener(alpha[t,p], tau[t,p], beta[t,p], delta[t,p])
    
    % params = alpha, tau, beta, delta, px, lr
   
%% Add paths for helper scripts and load in the master data struct and truncate it so there are only 3 subjects
environment = getenv('HOME');

% On Sherlock
if environment == '/home/users/mcsnyder'
    addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS'));
	addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/RL_DDM'));
	addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/mfit'));
    addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data'));
	base_path = '/oak/stanford/groups/menon/';

% On PigglyWig*2
elseif environment == '/home/scsnl'
	addpath(genpath('/home/scsnl/madeleinesnvzyder'));
	addpath(genpath('/home/scsnl/madeleinesnyder/RL_DDM'));
	addpath(genpath('/home/scsnl/madeleinesnyder/mfit'));
	addpath(genpath('/home/scsnl/madeleinesnyder/madeleine_data'));
	base_path = '/home/scsnl/madeleinesnyder/madeleine_data';
end 

addpath(genpath(strcat(base_path,'software/jags/4.3.0/lib/JAGS/modules-4')));
addpath(genpath(strcat(base_path,'software/jags')));
Asupp = STEP2_load_bart_supp(strcat(base_path,'projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data/master_qlearning_A.csv'));
A = STEP1_load_bart(strcat(base_path,'projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data/master_qlearning_A.csv'));

% Recode the reaction times to reflect choice and rt save as variable x
for i = 1:length(Asupp)
    for j = 1:length(Asupp(i).rt)
        if Asupp(i).c(j) == 1
            Asupp(i).x(j) = Asupp(i).rt(j)/1000;
        elseif Asupp(i).c(j) == 0
            Asupp(i).x(j) = -Asupp(i).rt(j)/1000;
        end
    end
end

% Makes subje index column 
for i=1:length(Asupp)
    Asupp(i).sidx = ones(length(Asupp(i).rt),1)*i;
end

% Recode HE
for i=1:length(Asupp)
    Asupp(i).HET = ones(length(Asupp(i).rt),1);
end

for j=1:length(Asupp)
    for i=1:length(Asupp(j).rt)
        if Asupp(j).HE(i) == '5'
            Asupp(j).HET(i) = '1';
        end 
        if Asupp(j).HE(i) == '3'
            Asupp(j).HET(i) = '2';
        end
        if Asupp(j).HE(i) == '0'
            Asupp(j).HET(i) = '2';
        end
    end 
end

for i=1:length(Asupp)
    AA(i).x = transpose(Asupp(i).x);
    AA(i).c = Asupp(i).c;
    AA(i).sidx = Asupp(i).sidx;
    AA(i).trial = Asupp(i).trial;
    AA(i).pumpidx = Asupp(i).s;
    AA(i).HE = Asupp(i).HET;
end

data = []; % preallocate space for the big data array

D = struct2cell(AA); % Makes AA1{1:370}, where AA1{1} is [69x1], AA1{2}
for j = 1:length(Asupp) % 1 to 74
    lens(j) = length(D{j*6});
end
slens = sum(lens); % get final lenght of big data vector 
ctr = 1 % initialize counter
T = zeros(slens,6); % preallocate space for the concatinated arrays
for j = 1:length(Asupp) % 1 to 74
    Atemp = vertcat(D{j*6-5:j*6}); % Create temporary array that holds the 6-entry swath of data starting from j*6 (index by 6)
    Atemp = reshape(Atemp,[lens(j),6]); % reshapes the [(length(D{j*5}*5) x 1] vector into a [(D{j*5}) x 5] array
    T(ctr:ctr+lens(j)-1,:) = Atemp;
    ctr = ctr + lens(j); % update the ctr so you know where to put the next array
end

datastruct = struct('x',T(:,1),'c',T(:,2),'sidx',T(:,3),'trial',T(:,4),'pumpidx',T(:,5),'lens',lens,'HE',T(:,6));

% MCMC parameters for JAGS
nchains=3; % How Many Chains?
nburnin=500; % How Many Burn-in Samples?
nsamples=1000; % How Many Recorded Samples?

% Initialize values all latent variables in all chains
clearvars init0
for i=1:nchains
    S.alpha = ones(74,2)*10; % Initial response threshold
    S.tau = ones(74,1)*0.001; % Initial non-decision time 
    S.beta = ones(74,1)*0.5; % Initial bias
    S.delta = ones(74,1)*0.1; % Initial drift rate
    S.px = T(:,1); % Initialze predicted values
    init0(i) = S; 
end

%%
% Use JAGS to sample
tic
doparallel = 0;
fprintf( 'Running JAGS\n' );
[samples, stats ] = matjags(datastruct,fullfile('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_scripts/wiener4.txt'),init0,'doparallel',doparallel,'nchains',nchains,'nburnin', nburnin,'nsamples', nsamples,'thin', 1,'dic', 1,'monitorparams', {'alpha','tau','beta','delta','px'},'savejagsoutput' , 1 ,'verbosity' , 1 ,'cleanup' , 0  );
toc





