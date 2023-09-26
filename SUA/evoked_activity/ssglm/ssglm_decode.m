%% Decode rewarded/non-rewarded trials
% Decode whether a trial was rewarded or not based on single unit activity
% and an estimated SS-GLM
%
% Compare to stationary model decoding
%
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
spike_ind = 1;
n_seg = 50;

session_info = setup(session);
save_folder = session_info.save_folder;

%% 
task_variant = 'rewarded';
variant_names = {'Incorrect','Correct'};
variants = [0,1]; 

%%
full_model_file = [save_folder,session,'_unit',num2str(spike_ind),'_',...
    task_variant,'_','R',num2str(n_seg)];
load(full_model_file)
stationary = load([full_model_file,'_stationary']);

%% Extract session information:

[session_info] = setup(session);
fs = session_info.fs_spike;
infusion_start = session_info.pharmacology.infusion_start;
infusion_end = session_info.pharmacology.infusion_end;
load('trial_info.mat')

%% Extract task information:

[spikes_all, time, task, trials, trial_times, correct_trials, ...
    complete_trials] = evoked_spike(session,spike_ind,task_variant);

trial_dur = time(trials(2,1))-time(trials(1,1));
trial_time = time(trials(1,1):trials(2,1))-time(trials(1,1));

l_seg = trial_dur/n_seg;
n_models = length(variants);
n_trials = length(trial_times);

%% Downsample spike data:

spikes_all_ms = zeros(round(trial_dur*1000),n_trials);

s_time = 0:1/fs:(length(spikes_all)-1)/fs;
t_ds = 0.001; % downsample spikes to ms 
s_ds = fs*t_ds;
n_s = round(trial_dur/t_ds);

time_ds = zeros(n_s,1);

for n = 1:n_s
    ind = ((n-1)*s_ds+1):(n*s_ds);
    time_ds(n) = mean(s_time(ind));
    spikes_all_ms(n,:) = sum(spikes_all(ind,:));
    spikes_all_ms(n,spikes_all_ms(n,:)>1) = 1;
end
Y = spikes_all_ms;

%% Remove outlier trials
% Noise can be improperly sorted into units. Remove trials corrupted by
% noise

for n = 1:2
    trial_outliers = remove_spike_outliers(Y,120);
    Y(:,trial_outliers) = [];
    trial_times(:,trial_outliers) = [];
    complete_trials(trial_outliers) = [];
    correct_trials(trial_outliers) = [];
    task(trial_outliers) = [];
    n_trials = size(Y,2);
end

%% Separate training and test trials

trials_keep = 1:n_trials; % Change if any trials are not to be included
est_trials = 1:2:n_trials; % Trials that were used for estimation

n_trials = length(trials_keep);
Y = Y(:,trials_keep);
n_Y = length(Y);
trial_times = trial_times(:,trials_keep);
complete_trials = complete_trials(trials_keep);
task = task(trials_keep);
complete_trials(est_trials) = 0; 
% ^ decoding will be run on complete trials

variant_trials = zeros(n_trials,n_models);
for n = 1:n_models
    variant_trials(:,n) = (task==variants(n)).*complete_trials;
end

n_trials_m = sum(variant_trials);
trial_mat = variant_trials==1; % logical version of variant_trials

%% Compute spike history + task basis
% This section uses the nSTAT toolbox

delta = time_ds(2)-time_ds(1);
fitType = 'poisson';
windowTimes = [0,1,5,10,15,20,30,40,50,100]/1000;

J = length(windowTimes)-1;
A = eye(n_seg,n_seg);

dN = Y';
K=size(dN,1);
minTime=0;
maxTime=(size(dN,2)-1)*delta;
histObj = History(windowTimes,minTime,maxTime);
Hk = cell(1,K);
spike_hist = zeros(n_s,n_trials,J);

for k=1:K
    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    nst{k}.setMinTime(minTime);
    nst{k}.setMaxTime(maxTime);
    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
    spike_hist(:,k,:) = Hk{k};
end


basisWidth = (maxTime-minTime)/n_seg;
sampleRate=1/delta;
unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,...
    minTime,maxTime,sampleRate);
basisMat = unitPulseBasis.data;
R = n_seg;


%% Decoding
% Use Monte Carlo to compute confidence intervals

trials_decode = complete_trials;
trial_list = find(trials_decode);
n_decode = sum(trials_decode);
models = [1,2];

dev = zeros(n_decode,2);
n_mc = 1000;
accuracy = zeros(n_decode,n_mc);
decoded = zeros(n_decode,n_mc);
xK_mc = cell(2,1);

dev_stat = zeros(n_decode,2);
accuracy_stat = zeros(n_decode,n_mc);
decoded_stat = zeros(n_decode,n_mc);
xK_mc_stat = cell(2,1);

correct_trials_decode = correct_trials(trials_decode==1);


ci_level = 0.05;


for b = 1:n_mc
    for n = 1:n_decode
        tic
        trial_ind = trial_list(n);

        mu_d = zeros(n_Y,n_models);
        mu_d_null = zeros(n_Y,n_models);

        for d=1:n_models
            m = models(d);
            model_ind = m;
            dN = Y(:,trial_ind)';
            Hk_trial=Hk{trial_ind};
            % Sample the distribution of xK:
            xK_mc{m} = mvnrnd(xK{m}(:,trial_ind),WK{m}(:,:,trial_ind),1)';
            % Compute history and stim effects:
            histEffect=exp(gammahat{m}*Hk_trial')';
            stimK=basisMat*xK_mc{m};
            stimEffect=exp(stimK);
            mu = histEffect.*stimEffect;
            % Compute the task encoding effect:
            stimEffect_scaled = stimEffect/mean(mu);
            % Scale the task encoding effect by the mean trial firing rate:
            trial_mu = sum(dN)/(n_s);
            stimEffect_scaled = stimEffect_scaled*trial_mu;
            % Compute the instantaneous firing rate during the trial:
            mu_d(:,d) = histEffect.*stimEffect_scaled;

            mu_adj = mu_d(:,d);
            y = dN';
            % Compute deviance:
            dev(n,d) = 2*sum(y.*(log((y+(y==0))./mu_adj))-(y-mu_adj));
            
            % Repeat the above steps for the stationary model:
            histEffect=exp(stationary.gammahat{m}*Hk_trial')';
            sigma_xK = stationary.stats{m}.se(1:n_seg);
            sigma_xK = diag(sigma_xK.^2);
            xK_mc_stat{m}= mvnrnd(stationary.xK{m},sigma_xK,1);
            stimK=basisMat*xK_mc_null{m}';
            stimEffect=exp(stimK);
            mu = histEffect.*stimEffect;
            stimEffect_scaled = stimEffect/mean(mu);
            trial_mu = sum(dN)/(n_s);
            stimEffect_scaled = stimEffect_scaled*trial_mu;
            mu_d_stat(:,d) = histEffect.*stimEffect_scaled;

            mu_adj = mu_d_stat(:,d);
            y = dN';
            dev_stat(n,d) = 2*sum(y.*(log((y+(y==0))./mu_adj))-(y-mu_adj));
        end
    end    
    
    % Compute decoding accuracy:
    dev_acc = zeros(size(dev,1),1);
    % Correct trials are accurately decoded if deviance for the incorrect
    % model is greater than deviance for the correct model
    dev_acc(correct_trials_decode==1) = ...
        dev(correct_trials_decode==1,1)...
        ./dev(correct_trials_decode==1,2);
    % Incorrect trials are accurately decoded if deviance for the correct
    % model is greater than deviance for the incorrect model
    dev_acc(correct_trials_decode==0) = ...
        dev(correct_trials_decode==0,2)...
        ./dev(correct_trials_decode==0,1);
    accuracy(:,b) = movmean(dev_acc>1,50); 
    
    dev_acc_stat = zeros(size(dev,1),1);
    dev_acc_stat(correct_trials_decode==1) = ...
        dev_stat(correct_trials_decode==1,1)...
        ./dev_stat(correct_trials_decode==1,2);
    dev_acc_stat(correct_trials_decode==0) = ...
        dev_stat(correct_trials_decode==0,2)...
        ./dev_stat(correct_trials_decode==0,1);
    [~,decoded_stat(:,b)] = min(dev_stat,[],2);
    accuracy_stat(:,b) = movmean(dev_acc_stat>1,50); 
end

%% Compute statistics

accuracy_est = sum(accuracy,2)/(n_mc);
accuracy_ci = zeros(size(accuracy_est,1),2);
for n = 1:n_decode
    accuracy_ci(n,1) = quantile(accuracy(n,:),ci_level/2);
    accuracy_ci(n,2) = quantile(accuracy(n,:),1-ci_level/2);
end

accuracy_est_stat = sum(accuracy_stat,2)/(n_mc);
accuracy_ci_stat = zeros(size(accuracy_est_stat,1),2);
for n = 1:n_decode
    accuracy_ci_stat(n,1) = quantile(accuracy_stat(n,:),ci_level/2);
    accuracy_ci_stat(n,2) = quantile(accuracy_stat(n,:),1-ci_level/2);
end

% Null decoding accuracy:
rewarded= correct_trials.*complete_trials;
p_rewarded = sum(rewarded)/sum(complete_trials);
null = p_rewarded.^2+(1-p_rewarded).^2;

% Compare periods:
before_trials = (trial_times(2,:)<infusion_start(1))';
after_trials = (trial_times(1,:)>infusion_end(end))';
inhibition_trials = (((trial_times(1,:)>infusion_end(end)).*...
    (trial_times(1,:)<(infusion_end(end)+10*60)))==1)';
recovery_trials = (((trial_times(1,:)>(infusion_end(end)+20*60)))==1)';

% Each period is characterized by trials in the middle of the period 
% (avoid edge effects)
before_decode= cumsum(trials_decode).*trials_decode.*before_tri%% Compute statistics

accuracy_est = sum(accuracy,2)/(n_mc);
accuracy_ci = zeros(size(accuracy_est,1),2);
for n = 1:n_decode
    accuracy_ci(n,1) = quantile(accuracy(n,:),ci_level/2);
    accuracy_ci(n,2) = quantile(accuracy(n,:),1-ci_level/2);
end

accuracy_est_stat = sum(accuracy_stat,2)/(n_mc);
accuracy_ci_stat = zeros(size(accuracy_est_stat,1),2);
for n = 1:n_decode
    accuracy_ci_stat(n,1) = quantile(accuracy_stat(n,:),ci_level/2);
    accuracy_ci_stat(n,2) = quantile(accuracy_stat(n,:),1-ci_level/2);
end

% Null decoding accuracy:
rewarded= correct_trials.*complete_trials;
p_rewarded = sum(rewarded)/sum(complete_trials);
null = p_rewarded.^2+(1-p_rewarded).^2;

% Compare periods:
before_trials = (trial_times(2,:)<infusion_start(1))';
after_trials = (trial_times(1,:)>infusion_end(end))';
inhibition_trials = (((trial_times(1,:)>infusion_end(end)).*...
    (trial_times(1,:)<(infusion_end(end)+10*60)))==1)';
recovery_trials = (((trial_times(1,:)>(infusion_end(end)+20*60)))==1)';

% Each period is characterized by trials in the middle of the period 
% (avoid edge effects)
before_decode= cumsum(trials_decode).*trials_decode.*before_trials;
before_decode(before_decode==0) = [];
before_decode = before_decode(1:end-15); 
after_decode= cumsum(trials_decode).*trials_decode.*after_trials;
after_decode(after_decode==0) = []; 
after_decode = after_decode(15:end);

inhib_decode= cumsum(trials_decode).*trials_decode.*inhibition_trials;
inhib_decode(inhib_decode==0) = [];
inhib_decode = inhib_decode(15:end-15); 

recov_decode= cumsum(trials_decode).*trials_decode.*recovery_trials;
recov_decode(recov_decode==0) = [];
recov_decode = recov_decode(15:end); 

n_before_decode = length(before_decode);
accuracy_before = accuracy(before_decode,:);

accuracy_before_est = quantile(mean(accuracy(before_decode,:)),0.5);
accuracy_before_ci = zeros(1,2);
accuracy_before_ci(1) = quantile(mean(accuracy(before_decode,:)),...
    ci_level/2);
accuracy_before_ci(2) = quantile(mean(accuracy(before_decode,:)),...
    1-ci_level/2);

accuracy_inhibition_est = quantile(mean(accuracy(inhib_decode,:)),0.5);
accuracy_inhibition_ci = zeros(1,2);
accuracy_inhibition_ci(1) = quantile(mean(accuracy(inhib_decode,:)),...
    ci_level/2);
accuracy_inhibition_ci(2) = quantile(mean(accuracy(inhib_decode,:)),...
    1-ci_level/2);

accuracy_recovery_est = quantile(mean(accuracy(recov_decode,:)),0.5);
accuracy_recovery_ci = zeros(1,2);
accuracy_recovery_ci(1) = quantile(mean(accuracy(recov_decode,:)),...
    ci_level/2);
accuracy_recovery_ci(2) = quantile(mean(accuracy(recov_decode,:)),...
    1-ci_level/2);

accuracy_sighigh = accuracy_ci(:,1)>null;
accuracy_sighigh_after = zeros(n_decode,1);
accuracy_sighigh_after(after_decode) = accuracy_sighigh(after_decode);
first_recovery = find(accuracy_sighigh_after,1,'first');
first_trial_recovery = trials_decode(first_recovery);

accuracy_belowstat = accuracy_ci(:,2) < accuracy_ci_stat(:,1);
accuracy_abovestat = accuracy_ci(:,1) > accuracy_ci_stat(:,2);
accuracy_equalorabovestat = ~accuracy_belowstat;
perc_equalorabove = sum(accuracy_equalorabovestat)/n_decode;
perc_above = sum(accuracy_abovestat)/n_decode;

%% Plot results

figure
trial_t = mean(trial_times(:,trials_decode==1))/60;
plotConfInterval_cont(trial_t,...
    accuracy_est',accuracy_ci(:,1)', accuracy_ci(:,2)')
hold on
plotConfInterval_cont(trial_t,...
    accuracy_est_stat',accuracy_ci_stat(:,1)', accuracy_ci_stat(:,2)',1)
plot(trial_t([first_recovery,first_recovery]),[0,1],'k')

plot([start_infuse/60,end_infuse/60],[1,1],'k')
plot([min(trial_t),max(trial_t)],ones(2,1)*null,'--k')
ylim([0,1.03])

xlim([0,max(trial_t)])als;
before_decode(before_decode==0) = [];
before_decode = before_decode(1:end-15); 
after_decode= cumsum(trials_decode).*trials_decode.*after_trials;
after_decode(after_decode==0) = []; 
after_decode = after_decode(15:end);

inhib_decode= cumsum(trials_decode).*trials_decode.*inhibition_trials;
inhib_decode(inhib_decode==0) = [];
inhib_decode = inhib_decode(15:end-15); 

recov_decode= cumsum(trials_decode).*trials_decode.*recovery_trials;
recov_decode(recov_decode==0) = [];
recov_decode = recov_decode(15:end); 

n_before_decode = length(before_decode);
accuracy_before = accuracy(before_decode,:);

accuracy_before_est = quantile(mean(accuracy(before_decode,:)),0.5);
accuracy_before_ci = zeros(1,2);
accuracy_before_ci(1) = quantile(mean(accuracy(before_decode,:)),...
    ci_level/2);
accuracy_before_ci(2) = quantile(mean(accuracy(before_decode,:)),...
    1-ci_level/2);

accuracy_inhibition_est = quantile(mean(accuracy(inhib_decode,:)),0.5);
accuracy_inhibition_ci = zeros(1,2);
accuracy_inhibition_ci(1) = quantile(mean(accuracy(inhib_decode,:)),...
    ci_level/2);
accuracy_inhibition_ci(2) = quantile(mean(accuracy(inhib_decode,:)),...
    1-ci_level/2);

accuracy_recovery_est = quantile(mean(accuracy(recov_decode,:)),0.5);
accuracy_recovery_ci = zeros(1,2);
accuracy_recovery_ci(1) = quantile(mean(accuracy(recov_decode,:)),...
    ci_level/2);
accuracy_recovery_ci(2) = quantile(mean(accuracy(recov_decode,:)),...
    1-ci_level/2);

accuracy_sighigh = accuracy_ci(:,1)>null;
accuracy_sighigh_after = zeros(n_decode,1);
accuracy_sighigh_after(after_decode) = accuracy_sighigh(after_decode);
first_recovery = find(accuracy_sighigh_after,1,'first');
first_trial_recovery = trials_decode(first_recovery);

accuracy_belowstat = accuracy_ci(:,2) < accuracy_ci_stat(:,1);
accuracy_abovestat = accuracy_ci(:,1) > accuracy_ci_stat(:,2);
accuracy_equalorabovestat = ~accuracy_belowstat;
perc_equalorabove = sum(accuracy_equalorabovestat)/n_decode;
perc_above = sum(accuracy_abovestat)/n_decode;

%% Plot results

figure
trial_t = mean(trial_times(:,trials_decode==1))/60;
plotConfInterval_cont(trial_t,...
    accuracy_est',accuracy_ci(:,1)', accuracy_ci(:,2)')
hold on
plotConfInterval_cont(trial_t,...
    accuracy_est_stat',accuracy_ci_stat(:,1)', accuracy_ci_stat(:,2)',1)
plot(trial_t([first_recovery,first_recovery]),[0,1],'k')

plot([start_infuse/60,end_infuse/60],[1,1],'k')
plot([min(trial_t),max(trial_t)],ones(2,1)*null,'--k')
ylim([0,1.03])

xlim([0,max(trial_t)])
