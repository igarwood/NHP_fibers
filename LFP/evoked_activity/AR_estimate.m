%% Estimate autoregressive model with task covariates for LFP data
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
task_variant = 'rewarded';
n_seg = 60;

session_info = setup(session);
save_folder = session_info.save_folder;

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
fs = session_info.fs_lfp;
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
test_trials = 2:2:n_trials; % Trials to be held out from estimation

n_trials = length(trials_keep);
Y = Y(:,trials_keep);
trial_times = trial_times(:,trials_keep);
complete_trials = complete_trials(trials_keep);
task = task(trials_keep);
complete_trials(test_trials) = 0; 
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

%% Plot estimation trials
figure
modulation_color = [0.65 0.9 0.9210];
for n = 1:n_variants
    subplot(1,n_variants,n)
    trials_plot = trial_mat(:,n);
    
    Y_plot = [Y(:,trials_plot)];
    
    imagesc(time(1:n_s),1:size(Y_plot,2),Y_plot')
    hold on
    
    for j = 1:length(infusion_start)
        modulation_trials = ...
            get_modulation_trials(trial_times(:,trials_plot),...
            infusion_start,infusion_end);
        x_fill = [time(1),time(n_s),time(n_s),time(1)];
        y_fill = [modulation_trials(1),modulation_trials(1),...
            modulation_trials(end),modulation_trials(end)];
        box = patch(x_fill,y_fill,modulation_color,'FaceAlpha',...
            0.7,'linestyle','none');
    end
    axis xy
    colormap('gray')
    Y_flat = reshape(Y_plot,1,[]);
    caxis([quantile(Y_flat,0.025),quantile(Y_flat,0.975)])
    xlabel('Time (s)')
    ylabel('Trial #')
    title([variant_names{n},' trials'])
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot LFP for matched trials from each variant
% Plot the same number of trials for each variant
% Plot all trials (not just trials to be estimated)
task_plot = task;
task(complete_trials==0)=-1;
if ~strcmp(task_variant,'all')
    plot_lfp_variants(Y, time(1:n_s), trial_times, ...
        infusion_start, infusion_end, variants, task)
end
for n = 1:n_variants
    subplot(1,n_variants,n)
    title([variant_names{n},' trials'])
    xlabel('Time (s)')
    ylabel('Trial #')
end

%% Estimate AR models with and without task encoding
for m = 1:n_variants
    model_ind = m;
    n_trials_m = sum(trial_mat(:,model_ind));  

    Y_m = Y(:,trial_mat(:,model_ind));
    Y_flat = reshape(Y_m,[],1);
    trial_segs_flat = repmat(basisMat,n_trials_m,1);

    X = trial_segs_flat;
    Mdl_with_encoding = arima(2,0,0);
    Mdl_without_encoding = arima(2,0,0);

    idxpre = 1:2;
    idxest = 3:numel(Y_flat);

%     
    EstMdl_with_encoding{m} = estimate(Mdl_with_encoding,Y_flat(idxest),'Y0',...
        Y_flat(idxpre),'X',X(idxest,:));
    [E,V,L] =  infer(EstMdl_with_encoding{m},Y_flat(idxest),'Y0',...
        Y_flat(idxpre),'X',X(idxest,:));
    E_with_encoding{m} = E;
    V_with_encoding{m} = V;
    LL_with_encoding(m) = L;
    EstMdl_without_encoding{m} = estimate(Mdl_without_encoding,Y_flat(idxest),...
        'Y0',Y_flat(idxpre));
    [E,V,L] = infer(EstMdl_without_encoding{m},Y_flat(idxest),...
        'Y0',Y_flat(idxpre));
    E_without_encoding{m} = E;
    V_without_encoding{m} = V;
    LL_without_encoding(m) = L;
end
save([save_folder,session,'_lfp_',task_variant,...
    '_','R',num2str(R)],'E_with_encoding', 'E_without_encoding', ...
    'EstMdl_with_encoding','EstMdl_without_encoding', ...
    'LL_with_encoding','LL_without_encoding','Mdl_with_encoding', ...
    'Mdl_without_encoding','V_with_encoding', 'V_without_encoding');

%% Generate new LFP samples from estimated models
Y_sim = zeros(size(Y));

for m = 1:n_models
    trials = find(trial_mat(:,m));

    Y_m = Y(:,trials);
    
    Y_flat = reshape(Y_m,[],1);
    Y_sim_flat = zeros(size(Y_flat));
    
    trial_segs_flat = repmat(basisMat,n_trials_m,1);
    X = trial_segs_flat;
    idxpre = 1:2;
    idxest = 3:numel(Y_flat);
    Y_sim_flat(idxpre) = Y_flat(idxpre);
    
    [Y_sim_flat(idxest)] = simulate(EstMdl_with_encoding{m},...
        length(Y_flat)-2,'Y0',Y_flat(idxpre),'X',X(idxest,:));
    Y_sim(:,trials) = reshape(Y_sim_flat,size(Y_m,1),size(Y_m,2));
end

%% Plot average simulated response
figure
for m = 1:n_models
    subplot(1,n_models,m)
    trials_plot = trial_mat(:,m);

    Y_plot = [Y_sim(:,trials_plot)];
    sim = plot((0:(n_s-1))/fs,mean(Y_plot,2),'k')
    hold on
    Y_plot = [Y(:,trials_plot)];
    rec = plot((0:(n_s-1))/fs,mean(Y_plot,2),'k--')
    
    yl = ylim;
    title([variant_names{m},' trials'])
    if strcmp(task_variant,'rewarded')
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1-variants(m));
    else
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1);
    end
    
    legend([sim(1),rec(1)],'Simulated LFP', 'Recorded LFP','Location',...
        'Northwest')
    xlabel('Time since cue onset (seconds)')
end