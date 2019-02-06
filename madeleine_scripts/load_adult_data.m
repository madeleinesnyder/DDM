function AA = load_adult_data()

% Add the paths
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/RL_DDM'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/mfit'));
addpath(genpath('/oak/stanford/groups/menon/projects/mcsnyder/2018_RLDDM_JAGS/madeleine_data'));
base_path = '/oak/stanford/groups/menon/';
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
end
