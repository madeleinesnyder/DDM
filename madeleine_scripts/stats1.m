% Stats yay

% Looooad the stats data
kstats = load('../madeleine_results/kid_stats_5.mat')
astats = load('../madeleine_results/adult_stats_5.mat')

% Looooad the samples
ksamples = load('../madeleine_results/kid_samples_5.mat')
asamples = load('../madeleine_results/adult_samples_5.mat')

%%
% Load the real data
AA = load_adult_data()
KK = load_kid_data()

%%
% JUNK (redo this) 
ctr = 1
prob_correct_pos = zeros(108,74);
prob_correct_neg = zeros(108,74);
for j=1:74
    temp_sub_real = AA(j).x';
    temp = asamples.samples.px(:,:,ctr:ctr+length(temp_sub_real)-1);  
    temp_sub_sim = reshape(temp,[3000,length(temp_sub_real)]);
    for i=1:length(temp_sub_real)
        if temp_sub_real(1,i) > 0
            prob_correct_pos(i,j) = sum(temp_sub_sim(:,i) < 0)/3000;
        if temp_sub_real(1,i) < 0
            prob_correct_neg(i,j) = sum(temp_sub_sim(:,i) > 0)/3000;
        end
        end
    end 
end  

%%
% Calculate accuracy of simualted data
% Returns a matrix of the accumulated absolute differences between the
% actual data and the simulated data, for each choice for each subject.
accuracy = zeros(108,74);
ctr = 1;
for j=1:74
    temp_sub_real = AA(j).x';
    temp = asamples.px(:,:,ctr:ctr+length(temp_sub_real)-1); 
    temp_sub_sim = reshape(temp,[3000,length(temp_sub_real)]);
    for i=1:length(temp_sub_real)
        tt = temp_sub_real(1,i);
        accuracy(i,j) = sum(abs(repmat(tt,length(temp_sub_sim),1) - temp_sub_sim(:,i))/tt);
        ctr = length(temp_sub_real);
    end 
end 
Total_score = sum(sum(accuracy))

%%
% Calculate the MSE for each sample of the 3000
% Subject one test
ctr = 1
data_real = AA(1).x'
data_sim = asamples.samples.px(:,:,ctr:ctr+length(data_real)-1);
data_sim = reshape(data_sim,3000,length(data_real));
ctr = ctr + length(data_real)-1;
MSEs = zeros(3000,length(data_real));

% MSE = 1/n (sum(true_value - est_value)^2)
% Example with the first value of the first subject
tt = repmat(AA(1).x(1),3000,1);
a = (tt-data_sim(:,1)).^2;
b = sum(a)/3000;
% So apparently the MSE of the first value of the first subject is 0.5

% Scale this shit up to calcaulating all MSEs for the first subject 
MSE_sub1 = zeros(length(data_real),1);
for i=1:length(data_real)
    tt_temp = repmat(AA(1).x(i),3000,1);
    a_temp = (tt_temp-data_sim(:,i)).^2;
    b_temp = sum(a_temp)/3000;
    MSE_sub1(i) = b_temp;
end 

% Scale this up to all subjects 
ctr = 1
MSE = zeros(108,74);
for j=1:74
    data_sim = asamples.samples.px(:,:,ctr:ctr+length(AA(j).x)-1);
    data_sim = reshape(data_sim,3000,length(AA(j).x));
    ctr = ctr + length(AA(j).x)-1;
    for i=1:length(AA(j).x)
        tt_temp = repmat(AA(j).x(i),3000,1);
        a_temp = (tt_temp-data_sim(:,i)).^2;
        b_temp = sum(a_temp)/3000;
        MSE(i,j) = b_temp;
    end 
end 
 
        
%%  
ctr = 1
MSE = zeros(6010,3000);
for j=1:74
    data_sim = asamples.samples.px(:,:,ctr:ctr+length(AA(j).x)-1);
    data_sim = reshape(data_sim,3000,length(AA(j).x));
    for i=1:length(AA(j).x)
        tt_temp = repmat(AA(j).x(i),3000,1);
        a_temp = (tt_temp-data_sim(:,i)).^2;
        b_temp = sum(a_temp);
        MSE(i,j) = b_temp;
    end
end
 
 














