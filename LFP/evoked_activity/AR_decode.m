%% Decode rewarded/non-rewarded trials
% Decode whether a trial was rewarded or not based on LFP activity
% and an estimated AR model
%
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
n_seg = 60;

session_info = setup(session);
save_folder = session_info.save_folder;

%% 
task_variant = 'rewarded';
variant_names = {'Incorrect','Correct'};
variants = [0,1]; 

%%
full_model_file = [save_folder,session,'_lfp_',...
    task_variant,'_','R',num2str(n_seg)];
load(full_model_file)

%% Extract session information:

[session_info] = setup(session);
fs = session_info.fs_spike;
infusion_start = session_info.pharmacology.infusion_start;
infusion_end = session_info.pharmacology.infusion_end;
load('trial_info.mat')

%% Extract task information:

[lfp_all, time, task, trials, trial_times, correct_trials, ...
    complete_trials] = evoked_LFP(session,task_variant);

trial_dur = trial_times(2,1)-trial_times(1,1);

l_seg = trial_dur/n_seg;
n_models = length(variants);
n_trials = length(trial_times);

Y = lfp_all;
n_s = size(Y,1);

%% Separate training and test trials
trials_keep = 1:n_trials; % Change if any trials are not to be included
est_trials = 1:2:n_trials; % Trials used for estimation

n_trials = length(trials_keep);
Y = Y(:,trials_keep);
trial_times = trial_times(:,trials_keep);
complete_trials = complete_trials(trials_keep);
correct_trials = correct_trials(trials_keep);
task = task(trials_keep);
complete_trials(est_trials) = 0; 
% ^ glm will be estimated from complete trials

%% Set up trial phase basis
minTime = 0;
maxTime = trial_dur-1/fs;
basisWidth = trial_dur/n_seg;
unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,...
    minTime,maxTime,fs);
basisMat = unitPulseBasis.data;
% Matlab's 'estimate' function includes a constant term; remove one column
% of basisMat to avoid singular design matrix
basisMat = basisMat(:,1:end-1); 
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

correct_trials_decode = correct_trials(trials_decode==1);

ci_level = 0.05;


for b = 1:n_mc
    for n = 1:n_decode
        tic
        trial_ind = trial_list(n);

        Y_n = Y(:,trial_ind);
        Y_flat = reshape(Y_n,[],1);
        %trial_segs_flat = repmat(basisMat,n_trials_m,1);
        X = basisMat;

        idxpre = 1:2;
        idxest = 3:numel(Y_flat);
        % Estimate AR params for test trial:
        EstMdl_test = estimate(Mdl_with_encoding,Y_flat(idxest),'Y0',...
            Y_flat(idxpre),'display','off');%,'X',X(idxest,:),'display','off');

        EstMdl_sample = EstMdl_with_encoding{1};
        EstMdl_sample.AR = EstMdl_test.AR; 

        for d=1:n_models
            m = models(d);
            model_ind = m;
            
            mdl_results = summarize(EstMdl_with_encoding{m});
            sn_sample = randn(mdl_results.NumEstimatedParameters,1);
            param_sample = mdl_results.Table{:,1}+...
                sn_sample.*mdl_results.Table{:,2};
            EstMdl_sample.Constant = param_sample(1);
            EstMdl_sample.Beta = param_sample(4:(n_seg+2))';
            EstMdl_sample.Variance = param_sample(end);
            
            
            [E,V,L] =  infer(EstMdl_sample,Y_flat(idxest),'Y0',...
                        Y_flat(idxpre),'X',X(idxest,:));
            dev(n,d) = sum(E.^2);
            
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
    accuracy(:,b) = movmean(dev_acc>1,30); 
    
end

%% Compute statistics

accuracy_est = sum(accuracy,2)/(n_mc);
accuracy_ci = zeros(size(accuracy_est,1),2);
for n = 1:n_decode
    accuracy_ci(n,1) = quantile(accuracy(n,:),ci_level/2);
    accuracy_ci(n,2) = quantile(accuracy(n,:),1-ci_level/2);
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



%% Plot results

figure
trial_t = mean(trial_times(:,trials_decode==1))/60;
plotConfInterval_cont(trial_t,...
    accuracy_est',accuracy_ci(:,1)', accuracy_ci(:,2)')
hold on
plot(trial_t([first_recovery,first_recovery]),[0,1],'k')

plot([infusion_start/60,infusion_end/60],[1,1],'k')
plot([min(trial_t),max(trial_t)],ones(2,1)*null,'--k')
ylim([0,1.03])

xlim([0,max(trial_t)])