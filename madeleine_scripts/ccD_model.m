addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/RL_DDM'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/mfit'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data'));
base_path = '/oak/stanford/groups/menon/';
addpath(genpath(strcat(base_path,'software/jags/4.3.0/lib/JAGS/modules-4')));
addpath(genpath(strcat(base_path,'software/jags')));

load('ccD.mat')

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