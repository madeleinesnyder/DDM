%%
% DDM taking into account pump-wise delta updating and belief states
% Runs belief1.txt JAGS model

% Things in the data struct
% response-coded RT*
% Cumulative number of bust trials up to the t_th trial*
% if the t_th trial was a bust trial
% if the t_th trial was a cash out
% shortest RT for the ith subject*
% max nuber of pumps possible in the task*
% number of subject*
% max number of trials*

% RTACC is the response coded RT
% cummBust(i,t) = Cumulative number of bust trials upto the t_th trial
% CI(i,t) = If the t_th trial was bust
% safe(i,t) = if the t_th trial was a cashout
% minT(i) = Shortest RT for i_th subject
% MX = maximum number of pumps possible in the task
% nI = number of subjects
% nT = maximum number of trials
   
%% Add paths for helper scripts and load in the master data struct and truncate it so there are only 3 subjects
environment = getenv('HOME');
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/RL_DDM'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/mfit'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data'));
base_path = '/oak/stanford/groups/menon/';
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

% Define new Asupp fields (subject index, history effect, cumulative busts,
% safe? risky? minimum RT (minT), maximum possible pumps (MX), number of
% subjects (nI), max nmber of trials (nT)
for i=1:length(Asupp)
    Asupp(i).sidx = ones(length(Asupp(i).rt),1)*i;
    Asupp(i).HET = ones(length(Asupp(i).rt),1);
    Asupp(i).cummBust = ones(length(Asupp(i).rt),1);
    Asupp(i).safe = ones(length(Asupp(i).rt),1);
    Asupp(i).CI = ones(length(Asupp(i).rt),1);
    Asupp(i).minT = ones(length(Asupp(i).rt),1);
    Asupp(i).MX = ones(length(Asupp(i).rt),1);
    Asupp(i).nI = ones(length(Asupp(i).rt),1)*74;
    Asupp(i).nT = ones(length(Asupp(i).rt),1);
end

% Recode HE
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

% cummBust
for j=1:length(Asupp)
    bust_ctr = 0;
    for i=1:length(Asupp(j).go)
        if Asupp(j).go(i) == 1
            bust_ctr = bust_ctr+1;
        end
    Asupp(j).cummBust(i) = bust_ctr;
    end 
end

% min rt (minT)
for j = 1:length(Asupp)
    min_rt = min(Asupp(j).rt)
    Asupp(j).minT(:) = repmat(min_rt,length(Asupp(j)),1);
end

% Max number of pumps in task
for j = 1:length(Asupp)
    max_pumps = max(Asupp(j).s);
    Asupp(j).MX(:) = repmat(max_pumps,length(Asupp(j)),1);
end 

% Max number of trials in the task (nT)
for j = 1:length(Asupp)
    max_trials = max(Asupp(j).exp_trial);
    Asupp(j).nT(:) = repmat(max_trials,length(Asupp(j)),1);
end 
        

for i=1:length(Asupp)
    AA(i).RTACC = transpose(Asupp(i).x);
    AA(i).c = Asupp(i).c;
    AA(i).sidx = Asupp(i).sidx;
    AA(i).trial = Asupp(i).exp_trial;
    AA(i).pumpidx = Asupp(i).s;
    AA(i).HE = Asupp(i).HET;
    AA(i).cummBust = Asupp(i).cummBust;
    AA(i).minT = Asupp(i).minT;
    AA(i).MX = Asupp(i).MX;
    AA(i).nI = Asupp(i).nI;
    AA(i).nT = Asupp(i).nT;
end

%%
RTACC = zeros(74,35,12)*nan;
for i = 1:length(Asupp)
    ctr = 1;
    j = 1;
    for j = 1:max(Asupp(i).exp_trial)
        k = 1;
        for k = 1:sum(Asupp(i).exp_trial == j)
            RTACC(i,j,k) = Asupp(i).x(ctr);
            ctr = ctr+1;
        end
    end
end

%%
D = struct2cell(AA); % Makes AA1{1:370}, where AA1{1} is [69x1], AA1{2}
for j = 1:length(Asupp) % 1 to 74
    lens(j) = length(D{j*11});
end
slens = sum(lens); % get final lenght of big data vector 
ctr = 1 % initialize counter
T = zeros(slens,11); % preallocate space for the concatinated arrays
for j = 1:length(Asupp) % 1 to 74
    Atemp = vertcat(D{(j*11)-10:j*11}); % Create temporary array that holds the 6-entry swath of data starting from j*6 (index by 6)
    Atemp = reshape(Atemp,[lens(j),11]); % reshapes the [(length(D{j*5}*5) x 1] vector into a [(D{j*5}) x 5] array
    T(ctr:ctr+lens(j)-1,:) = Atemp;
    ctr = ctr + lens(j); % update the ctr so you know where to put the next array
end

ccD = struct('RTACC',RTACC)%,'safe',T(:,2),'sidx',T(:,3),'trial',T(:,4),'pumpidx',T(:,5),'lens',lens,'HE',T(:,6),'cummBust',T(:,7),'minT',T(:,8),'MX',12,'nI',74,'nT',35);

% Eliminate the bad rts
ccD.minT(find(ccD.minT < 0)) = nan;
ccD.minT = min(ccD.minT);

% MCMC parameters for JAGS
nChains=3; % How Many Chains?
nburnin=500; % How Many Burn-in Samples?
nsamples=1000; % How Many Recorded Samples?

% Initial Values
for i=1:nChains
     S.SS=0;
     S.RTACC=nan(ccD.nI,ccD.nT,ccD.MX);
     S.RTACC(isnan(ccD.RTACC))=0.5.*(1-2.*ccD.safe(isnan(ccD.RTACC)));
     S.pRTACC=S.RTACC;
     S.pBelief = 0.5.*ones(ccD.nI,ccD.nT);
     S.beta = 0.5.*ones(1,ccD.nI);
     S.delta0 = 3.*ones(1,ccD.nI);
     S.dA = zeros(ccD.nI,2);
     S.r = 0.5.*ones(1,ccD.nI);
     S.s = 2.*ones(1,ccD.nI);
     S.k = 3.*ones(ccD.nI,2);
     S.tau = 0.01.*ones(1,ccD.nI);
     S.gamma = 2.*ones(1,ccD.nI);
    %S.pRTACC=ccD.RTACC;
    %S.pRTACC(isnan(ccD.RTACC))=0.5.*(1-2.*ccD.safe(isnan(ccD.RTACC)));
    init0(i) = S;
end

%%
% Use JAGS to sample
tic
doparallel = 0;
fprintf( 'Running JAGS\n' );
[samples, stats ] = matjags(ccD,fullfile('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_scripts/belief1.txt'),init0,'doparallel',doparallel,'nchains',nChains,'nburnin', nburnin,'nsamples', nsamples,'thin', 1,'dic', 1,'monitorparams', {'dA','r','s','K','t0','k','beta','delta0','alpha','delta','tau','muT','muG','precT','precG','pBelief','gamma','omega','pRTACC'},'savejagsoutput' , 1 ,'verbosity' , 1 ,'cleanup' , 0  );
toc





