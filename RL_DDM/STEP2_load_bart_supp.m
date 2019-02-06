function data = STEP2_load_bart_supp(csvfile)
    
    % Load data from BART task.
    %
    % USAGE: data = load_bart_supp(csvfile i.e. '/home/scsnl/madeleinesnyder/madeleine_data/master_qlearning_A.csv')
    %
    % OUTPUTS:
    %   data - [S x 1] structure, where S is the number of subjects, with the following fields:
    %           .c - [N x 1] choices
    %           .r - [N x 1] rewards
    %           .s - [N x 1] states
    %           .isblow - [N x 1] explode trial indicator (1=Explode, 0=Not Explode)
    %           .rt - [N x 1] response times
    %		.trial - [N x 1] trial number
    %       .exp_trial - [N x 1] trial number CONDITIONED ON EXPERIMENTAL
    %       TRIALS ONLY
    %		.HE - [N x 1] History effect (5 = post-explode, 3 = post-cashout
    %		.PCPP - [N x 1] - precashout or prepump (prepump = 1 precashout = 2)
    %           .C - number of choice options
    %           .N - number of trials
    
    D = csvread(csvfile,1);
    
    % Remove the rows that have 0 for the reaction time 
    indices = find(D(:,6)==0);
    D(indices,:) = [];
    
    subs = unique(D(:,1));      % subjects
    
    for i = 1:length(subs)
        ix = D(:,1)==subs(i);
        data(i).go = D(ix,2);
        data(i).s = D(ix,3);
        data(i).c = D(ix,4);
        data(i).r = D(ix,5);
        data(i).rt = D(ix,6);
        data(i).trial = D(ix,7);
        data(i).HE = D(ix,8);
        data(i).PCPP = D(ix,9);
        data(i).C = 2;
        data(i).N = length(data(i).c);
        data(i).exp_trial = ones(6010,1);
    end
    
    % Make the experiemntal-only trials
    for i = 1:length(subs)
        ix = D(:,1)==subs(i);
        TT = D(ix,7); 
        if TT(1) ~= 1
            data(i).exp_trial = TT-TT(1)+1
        elseif TT(1) == 1
            data(i).exp_trial = data(i).trial;
        end
    end
end

