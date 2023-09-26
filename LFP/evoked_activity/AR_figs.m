%% Plot LFP observations and average response + confidence intervals 
% Prerequisite: estimated AR model for the corresponding task variant
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
spike_ind = 1;
task_variant = 'rewarded';
n_seg = 60;

session_info = setup(session);
save_folder = session_info.save_folder;
%%
load([save_folder,session,'_lfp_',task_variant,...
    '_','R',num2str(R)]);
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

%% Separate trials into variants
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

%% Monte carlo simulation of LFP
% Note that storing all samples + parallelizing takes a lot of memory
% Could be optimized to reduce memory burden
n_mc = 1000;
Y_sim = zeros(size(Y,1),size(Y,2),n_bs);
parfor n = 1:n_mc
    tic
    Y_sim_loop = zeros(size(Y));
    for m = 1:n_models
        trials = find(trial_mat(:,m));
        n_trials_m = length(trials);  

        seed_sample = randsample(trials,n_trials_m);
        Y_sample = Y(:,seed_sample);
        
        Y_flat = reshape(Y_sample,[],1);
        Y_sim_flat = zeros(size(Y_flat));

        trial_segs_flat = repmat(basisMat,n_trials_m,1);
        X = trial_segs_flat;
        idxpre = 1:2;
        idxest = 3:numel(Y_flat);
        Y_sim_flat(idxpre) = Y_flat(idxpre);

        [Y_sim_flat(idxest)] = simulate(EstMdl_with_encoding{m},...
            length(Y_flat)-2,'Y0',Y_flat(idxpre),'X',X(idxest,:));
        Y_sim_loop(:,trials) = reshape(Y_sim_flat,n_s,n_trials_m);
    end
    Y_sim(:,:,n) = Y_sim_loop;
    toc
end

%% Extract quantiles
lfp_est = zeros(n_s,n_variants);
lfp_cf = zeros(n_s,2,n_variants);
lfp_true = zeros(n_s,n_variants);
for l = 1:n_models
    trials = find(trial_mat(:,l));
    
    for n = 1:n_s
%         lfp_est(n,l) = median(Y_sim(n,:,l));
        mean_sim = squeeze(mean(Y_sim(n,trials,:),2));
        lfp_cf(n,1,l) = quantile(mean_sim,0.025);
        lfp_cf(n,2,l) = quantile(mean_sim,0.975);
        lfp_true(n,l) = sum(Y(n,trials))/length(trials);
    end

end

%%
figure
perc_within_bounds = zeros(1,n_variants);
ax = [];
for n = 1:n_variants
    subplot(1,n_variants,n)
    plotConfInterval_cont(time(1:n_s),...
        lfp_true(:,n)',lfp_cf(:,1,n)', lfp_cf(:,2,n)')
    hold on
    
    perc_within_bounds(n) = sum((lfp_true(:,l)<lfp_cf(:,2,n)).*...
        (lfp_true(:,n)<lfp_cf(:,2,n)))./(n_s);
    ax = [ax,gca]
    if strcmp(task_variant,'rewarded')
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1-variants(n));
    else
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1);
    end
end
linkaxes(ax,'y');