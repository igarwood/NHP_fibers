%% Goodness of fit assessment for LFP fit to AR models with task covariants
% Test the GOF for task coefficients 
% AR parameters are reestimated for each trial
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

variant_trials = zeros(n_trials,n_models);
for n = 1:n_models
    variant_trials(:,n) = (task==variants(n)).*complete_trials;
end

n_trials_m = sum(variant_trials);
trial_mat = variant_trials==1; % logical version of variant_trials


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


%% LFP inference from estimated models
% Use Monte Carlo to compute confidence intervals

E_inf = cell(1,n_models);


for m = 1:n_models
    model_ind = m;
    n_trials_m = sum(trial_mat(:,model_ind));  

    Y_m = Y(:,trial_mat(:,model_ind));
    for n = 1:n_trials_m
        Y_n = Y_m(:,n);
        Y_flat = reshape(Y_n,[],1);
        X = basisMat;

        idxpre = 1:2;
        idxest = 3:numel(Y_flat);
        % Estimate AR params for test trial:
        EstMdl_test_estAR = estimate(Mdl_with_encoding,Y_flat(idxest),'Y0',...
            Y_flat(idxpre),'display','off');

        EstMdl_test = EstMdl_with_encoding{m};
        EstMdl_test.AR = EstMdl_test_estAR.AR; 
        
        [E,V,L] =  infer(EstMdl_test,Y_flat(idxest),'Y0',...
                    Y_flat(idxpre),'X',X(idxest,:));
        ind_a = (n-1)*(n_s-2)+1;
        ind_b = n*(n_s-2);
        E_inf{m}(ind_a:ind_b) = E;
    end
end


%% Plot ACF
% Number of lags correspond to ACF computed across 10 trials

figure
for m=1:n_models
    subplot(1,n_models,m)
    residual = E_inf{m};
    residual_var = var(E_inf{m});
    residual_standard = residual/sqrt(residual_var);
    [test_acf,test_acf_lags] = autocorr(norminv(residual_standard),...
        'NumLags',60000);
    
    rand_sample = randsample(2:60000,1000);
    plot(test_acf_lags(rand_sample),test_acf(rand_sample),'.')
    hold on
    plot([0,60000],zeros(2,1),'k')
    % CI are for gaussian white noise:
    plot([0,60000],ones(2,1)*2*sqrt(1/length(residual_standard)),'k--')
    plot([0,60000],ones(2,1)*-2*sqrt(1/length(residual_standard)),'k--')
    ylim([-0.01,0.03])
end