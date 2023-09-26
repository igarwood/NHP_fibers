%% Goodness of fit assessment for ssglm models
% Prerequisite: for a given trial variant, estimate SSGLM models 
% with history terms for R = 1 and R = 50
% without history terms for R = 1 and R = 50
% stationary (constant parameters across trials)
%
% This script is computationally intensive and takes a while to run
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
spike_ind = 1;
task_variant = 'rewarded';
n_seg = 50;
n_seg_null = 1;

session_info = setup(session);
save_folder = session_info.save_folder;

model_root = [save_folder,session,'_unit',num2str(spike_ind),'_',...
    task_variant,'_R'];
load([model_root,num2str(n_seg)]);
null = load([model_root,'1_nohist']); % null model has no history terms and does not vary within trials
no_encoding = load([model_root,'1']); % no_encoding model assumes rate is constant within trials
no_hist = load([model_root,num2str(n_seg),'_nohist']); % no_hist model does not include spike history terms
stationary = load([model_root,num2str(n_seg),'_stationary']); % stationary model does not vary across trials
%%
if strcmp(task_variant,'rewarded')    
    variant_names = {'Incorrect','Correct'};
    variants = [0,1]; 
elseif strcmp(task_variant,'sample')
    variant_names = {'Sample a','Sample b','Sample c'};
    variants = [1,2,3];
elseif strcmp(task_variant,'choice')
    variant_names = {'Choice a','Choice b','Choice c'};
    variants = [1,2,3];
elseif strcmp(task_variant,'lr')
    variant_names = {'Left','Right'};
    variants = [0,1];
elseif strcmp(task_variant,'all')
    variant_names = {'All'};
    variants = [0];
end
n_variants = length(variants);

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

%% Specify trials belonging to each task variant

trials_keep = 1:n_trials; % Change if any trials are not to be included


n_trials = length(trials_keep);
Y = Y(:,trials_keep);
trial_times = trial_times(:,trials_keep);
complete_trials = complete_trials(trials_keep);
task = task(trials_keep);

variant_trials = zeros(n_trials,n_models);
for n = 1:n_models
    variant_trials(:,n) = (task==variants(n)).*complete_trials;
end

n_trials_m = sum(variant_trials);
trial_mat = variant_trials==1; % logical version of variant_trials

%% Compute spike history
% This section uses the nSTAT toolbox

delta = time_ds(2)-time_ds(1);
fitType = 'poisson';
windowTimes = [0,1,5,10,15,20,30,40,50,100]/1000;

J = length(windowTimes)-1;
A = eye(n_seg,n_seg);
A_null=eye(1,1);

dN = Y';
K=size(dN,1);
minTime=0;
maxTime=(size(dN,2)-1)*delta;
histObj = History(windowTimes,minTime,maxTime);
Hk = cell(1,K);
Hk_null = cell(1,K);
spike_hist = zeros(n_s,n_trials,J);

for k=1:K
    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    nst{k}.setMinTime(minTime);
    nst{k}.setMaxTime(maxTime);
    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
    spike_hist(:,k,:) = Hk{k};
    Hk_null{k} = zeros(size(Y_m,2),1);
end


basisWidth = (maxTime-minTime)/n_seg;
sampleRate=1/delta;
unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,...
    minTime,maxTime,sampleRate);
basisMat = unitPulseBasis.data;
R = n_seg;
R_null = n_seg_null;
%% Data structures
xK_inf = cell(n_models,1);
WK_inf = cell(n_models,1);
Wku_inf = cell(n_models,1);
logll_inf = cell(n_models,1);
fitResults_inf = cell(n_models,1);
fR_inf = cell(n_models,1);
stimulus_inf = cell(n_models,1);
stimCIs_inf = cell(n_models,1);


xK_null = cell(n_models,1);
WK_null = cell(n_models,1);
Wku_null = cell(n_models,1);
logll_null = cell(n_models,1);
fitResults_null = cell(n_models,1);
fR_null = cell(n_models,1);

xK_no_encoding = cell(n_models,1);
WK_no_encoding = cell(n_models,1);
Wku_no_encoding = cell(n_models,1);
logll_no_encoding = cell(n_models,1);
fitResults_no_encoding = cell(n_models,1);
fR_no_encoding = cell(n_models,1);

xK_no_hist = cell(n_models,1);
WK_no_hist = cell(n_models,1);
Wku_no_hist = cell(n_models,1);
logll_no_hist = cell(n_models,1);
fitResults_no_hist = cell(n_models,1);
fR_no_hist = cell(n_models,1);

%% For all trials and ssglms, estimate rate and compute fit statistics

for m = 1:n_models
    model_ind = m;
    Y_m = Y(:,trial_mat(:,model_ind))';
    K=size(Y,2);
    
    spikes_task = Y_m;
    dN = zeros(length(complete_trials),size(spikes_all_ms,1));
    dN(trial_mat(:,model_ind),:) = spikes_task;
    missing_trials = ~trial_mat(:,model_ind);
    
    % Estimate rate + other ssglm parameters ------------------------------
    xK0=xK{m}(:,1);

    [xK_inf{m},WK_inf{m},Wku_inf{m},logll_inf{m},SumXkTerms,sumPPll]= ...
        DecodingAlgorithms_IG.PPSS_EStep_missing(A,Qhat{m},xK0,dN,...
        missing_trials,Hk,fitType,delta,gammahat{m},n_seg);

    xK0=xK_inf{m}(:,1);
    [xK_inf{m},WK_inf{m},Wku_inf{m},logll_inf{m},SumXkTerms,sumPPll]= ...
        DecodingAlgorithms_IG.PPSS_EStep_missing(A,Qhat{m},xK0,dN,...
        missing_trials,Hk,fitType,delta,gammahat{m},n_seg);
    xK0=xK_inf{m}(:,1);
    
    xK0_null=null.xK{m}(:,1);
    [xK_null{m},WK_null{m},Wku_null{m},logll_null{m},SumXkTerms_null,...
        sumPPll_null]= DecodingAlgorithms_IG.PPSS_EStep_missing(A_null,...
        null.Qhat{m},xK0_null,dN,missing_trials,Hk_null,fitType,delta,...
        null.gammahat{m},n_seg_null);
    xK0_null=xK_null{m}(:,1);
    [xK_null{m},WK_null{m},Wku_null{m},logll_null{m},SumXkTerms_null,...
        sumPPll_null]= ...
        DecodingAlgorithms_IG.PPSS_EStep_missing(A_null,null.Qhat{m},...
        xK0_null,dN,missing_trials,Hk_null,fitType,delta,...
        null.gammahat{m},n_seg_null);
    xK0_null=xK_null{m}(:,1);
    
    xK0_no_encoding=no_encoding.xK{m}(:,1);
    [xK_no_encoding{m},WK_no_encoding{m},Wku_no_encoding{m},...
        logll_no_encoding{m},SumXkTerms_no_encoding,...
        sumPPll_no_encoding]= DecodingAlgorithms_IG.PPSS_EStep_missing...
        (A_null,no_encoding.Qhat{m},xK0_no_encoding,dN,missing_trials,...
        Hk,fitType,delta,no_encoding.gammahat{m},n_seg_null);
    xK0_no_encoding=xK_no_encoding{m}(:,1);
    [xK_no_encoding{m},WK_no_encoding{m},Wku_no_encoding{m},...
        logll_no_encoding{m},SumXkTerms_no_encoding,...
        sumPPll_no_encoding]= DecodingAlgorithms_IG.PPSS_EStep_missing...
        (A_null,no_encoding.Qhat{m},xK0_no_encoding,dN,missing_trials,...
        Hk,fitType,delta,no_encoding.gammahat{m},n_seg_null);
    xK0_no_encoding=xK_no_encoding{m}(:,1);
    
    xK0_no_hist=no_hist.xK{m}(:,1);
    [xK_no_hist{m},WK_no_hist{m},Wku_no_hist{m},...
        logll_no_hist{m},SumXkTerms_no_hist,...
        sumPPll_no_hist]= DecodingAlgorithms_IG.PPSS_EStep_missing...
        (A,no_hist.Qhat{m},xK0_no_hist,dN,missing_trials,...
        Hk_null,fitType,delta,no_hist.gammahat{m},n_seg);
    xK0_no_hist=xK_no_hist{m}(:,1);
    [xK_no_hist{m},WK_no_hist{m},Wku_no_hist{m},...
        logll_no_hist{m},SumXkTerms_no_hist,...
        sumPPll_no_hist]= DecodingAlgorithms_IG.PPSS_EStep_missing...
        (A,no_hist.Qhat{m},xK0_no_hist,dN,missing_trials,...
        Hk_null,fitType,delta,no_hist.gammahat{m},n_seg);
    xK0_no_hist=xK_no_hist{m}(:,1);

    % Compute fit statistics ----------------------------------------------
    tic
    logllobs = logll_inf{m}+R*K*log(2*pi)+K/2*log(det(diag(Qhat{m})))+ ...
        1/2*trace(diag(Qhat{m})\SumXkTerms);
    fitResults_inf{m} = DecodingAlgorithms_IG.prepareEMResults_missing...
        (fitType,1,dN,Hk,xK_inf{m},WK_inf{m},Qhat{m},gammahat{m},...
        windowTimes,delta,[],logllobs);
    fR_inf{m} = fitResults_inf{m}.toStructure;
    toc
    
    tic
    logllobs_no_encoding = logll_no_encoding{m}+R_null*K*log(2*pi)+...
        K/2*log(det(diag(no_encoding.Qhat{m})))+ ...
        1/2*trace(diag(no_encoding.Qhat{m})\SumXkTerms_no_encoding);
    fitResults_no_encoding{m} = ...
        DecodingAlgorithms_IG.prepareEMResults_missing...
        (fitType,1,dN,Hk,xK_no_encoding{m},WK_no_encoding{m},...
        no_encoding.Qhat{m},no_encoding.gammahat{m},windowTimes,delta,...
        [],logllobs_no_encoding);
    fR_no_encoding{m} = fitResults_no_encoding{m}.toStructure;
    toc
    
    tic
    logllobs_no_hist = logll_no_hist{m}+R*K*log(2*pi)+K/2*...
        log(det(diag(no_hist.Qhat{m})))+ ...
        1/2*trace(diag(no_hist.Qhat{m})\SumXkTerms_no_hist);
    fitResults_no_hist{m} = ...
        DecodingAlgorithms_IG.prepareEMResults_missing(fitType,1,dN,...
        Hk_null,xK_no_hist{m},WK_no_hist{m},no_hist.Qhat{m},...
        no_hist.gammahat{m},[],delta,[],logllobs_no_hist);
    fR_no_hist{m} = fitResults_no_hist{m}.toStructure;
    toc
    
    tic 
    logllobs_null = logll_null{m}+R_null*K*log(2*pi)+K/2*...
        log(det(diag(null.Qhat{m})))+ ...
        1/2*trace(diag(null.Qhat{m})\SumXkTerms_null);
    fitResults_null{m} = ...
        DecodingAlgorithms_IG.prepareEMResults_missing(fitType,1,dN,...
        Hk_null,xK_null{m},WK_null{m},null.Qhat{m},null.gammahat{m},[],...
        delta,[],logllobs_null);
    fR_null{m} = fitResults_null{m}.toStructure;
    toc
    
end
%% For all trials, compute rate and fit statistics for the stationary model

for m = 1:n_models
    model_ind = m;
    n_trials_m = sum(trial_mat(:,model_ind));
    
    init_inds = find(trial_mat(:,model_ind));
    Y_m = Y(:,init_inds);

    Y_flat = reshape(Y(:,init_inds),[],1);
    spike_hist_flat = reshape(spike_hist(:,init_inds,:),[],J);
    trial_segs_flat = repmat(basisMat,length(init_inds),1);

    X_flat = [trial_segs_flat,spike_hist_flat];
    X = reshape(X_flat,[],n_trials_m,J+R);
    
    % Estimate rate + other glm parameters --------------------------------
    B_stationary = [stationary.xK{m};stationary.gammahat{m}'];
    n_Y = size(Y,1);
    [est_lamb_trials,est_lamb_trials_ci] = estimate_lambda(X,...
        B_stationary,stationary.stats{model_ind},n_Y);
    xK_stat{m} = repmat(B_stationary(1:n_seg),1,n_trials_m);
    gammahat_stat{m} = B_stationary(n_seg+1:end);
    
    % Compute fit statistics ----------------------------------------------
    Z = [];
    U = [];
    for n = 1:n_trials_m
        [z,u] = time_rescaling_theorem(est_lamb_trials(:,n),Y_m(:,n));
        Z = [Z;z];
        U = [U;u];
    end
    fR_stat{m}.Z = Z;
    fR_stat{m}.U = U;
    [~, pValue, ~, ~] = kstest(norminv(U));
    fR_stat{m}.KSStats.pValue = pValue;
end
%% KS-plots

figure
model_labels = {'Rewarded trials',...
    'Non-rewarded trials'};
for m = 1:2
    subplot(1,2,m)
    [f,x] = ecdf(fR_inf{m}.U);
    plot(x,f)  
    hold on
    
    [f,x] = ecdf(fR_null{m}.U);
    plot(x,f);  
    [f,x] = ecdf(fR_no_encoding{m}.U);
    plot(x,f);  
    [f,x] = ecdf(fR_no_hist{m}.U);
    plot(x,f); 
    [f,x] = ecdf(fR_stat{m}.U);
    plot(x,f);  
    n_int = length(fR_inf{m}.KSStats.KSSorted);
%     plot(fR{m}.KSStats.xAxis,fR{m}.KSStats.KSSorted)
    plot(fR_inf{m}.KSStats.xAxis,...
        fR_inf{m}.KSStats.xAxis+1.36/sqrt(n_int),'k--');
    plot(fR_inf{m}.KSStats.xAxis,...
        fR_inf{m}.KSStats.xAxis-1.36/sqrt(n_int),'k--');
    title(model_labels{m})
%     legend(['p-value = ', num2str(fR_inf{m}.KSStats.pValue)],...
%         'Location','southeast');
    %['p-value = ', num2str(fR{m}.KSStats.pValue)],...
    xlim([0,1])
    ylim([0,1])
    ylabel('Empirical CDF')
    xlabel('Model CDF')
end
legend({'full model', 'null model', 'no encoding', 'no history',...
    'stationary'})


%% ACF plot

figure
model_labels = {'Rewarded trials',...
    'Non-rewarded trials'};
for m = 1:2
    subplot(1,2,m)
    u = norminv(fR_inf{m}.U);
    u(isinf(u)) = [];
    [test_acf,test_acf_lags] = autocorr(u,'NumLags',1000);
    plot(test_acf_lags(2:end),test_acf(2:end),'.')
    hold on
    title(model_labels{m})
    %ylim([-0.1,0.1])
    xlabel('Lags')
    ylabel('ACF')
    
    u = norminv(fR_null{m}.U);
    u(isinf(u)) = [];
    [test_acf,test_acf_lags] = autocorr(u,'NumLags',1000);
    plot(test_acf_lags(2:end),test_acf(2:end),'.')
     
    u = norminv(fR_no_encoding{m}.U);
    u(isinf(u)) = [];
    [test_acf,test_acf_lags] = autocorr(u,'NumLags',1000);
    plot(test_acf_lags(2:end),test_acf(2:end),'.')
     
    u = norminv(fR_no_hist{m}.U);
    u(isinf(u)) = [];
    [test_acf,test_acf_lags] = autocorr(u,'NumLags',1000);
    plot(test_acf_lags(2:end),test_acf(2:end),'.')
    
    u = norminv(fR_stat{m}.U);
    u(isinf(u)) = [];
    [test_acf,test_acf_lags] = autocorr(u,'NumLags',1000);
    plot(test_acf_lags(2:end),test_acf(2:end),'.')
    
    plot([0,1000],zeros(2,1),'k')
    
    u = norminv(fR_inf{m}.U);
    % CI are for gaussian white noise:
    plot([0,1000],zeros(2,1),'k')
    plot([0,1000],ones(2,1)*2*sqrt(1/length(u)),'k--')
    plot([0,1000],ones(2,1)*-2*sqrt(1/length(u)),'k--')
    
end