function KK = load_kid_data()

addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/RL_DDM'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/mfit'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data'));
base_path = '/oak/stanford/groups/menon/';

Ksupp = STEP2_load_bart_supp(strcat(base_path,'projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data/master_qlearning_K.csv'));
K = STEP1_load_bart(strcat(base_path,'projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data/master_qlearning_K.csv'));

% Recode the reaction times to reflect choice and rt save as variable x
for i = 1:length(Ksupp)
    for j = 1:length(Ksupp(i).rt)
        if Ksupp(i).c(j) == 1
            Ksupp(i).x(j) = Ksupp(i).rt(j)/1000;
        elseif Ksupp(i).c(j) == 0
            Ksupp(i).x(j) = -Ksupp(i).rt(j)/1000;
        end
    end
end

% Makes subje index column 
for i=1:length(Ksupp)
    Ksupp(i).sidx = ones(length(Ksupp(i).rt),1)*i;
end

% Recode HE
for i=1:length(Ksupp)
    Ksupp(i).HET = ones(length(Ksupp(i).rt),1);
end

for j=1:length(Ksupp)
    for i=1:length(Ksupp(j).rt)
        if Ksupp(j).HE(i) == '5'
            Ksupp(j).HET(i) = '1';
        end 
        if Ksupp(j).HE(i) == '3'
            Ksupp(j).HET(i) = '2';
        end
        if Ksupp(j).HE(i) == '0'
            Ksupp(j).HET(i) = '2';
        end
    end 
end

for i=1:length(Ksupp)
    KK(i).x = transpose(Ksupp(i).x);
    KK(i).c = Ksupp(i).c;
    KK(i).sidx = Ksupp(i).sidx;
    KK(i).trial = Ksupp(i).trial;
    KK(i).pumpidx = Ksupp(i).s;
    KK(i).HE = Ksupp(i).HET;
end
end 
