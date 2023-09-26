%% Load spike data + SSGLM and create figures
% Prerequisite: run ssglm_estimate
% Fig 1: Raster plots of each trial type
% Fig 2: SSGLM rate estimation 
% Fig 3: Task encoding effect
%
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
spike_ind = 1;
task_variant = 'rewarded';
n_seg = 50;

session_info = setup(session);
save_folder = session_info.save_folder;

full_model_file = [save_folder,session,'_unit',num2str(spike_ind),'_',...
    task_variant,'_','R',num2str(n_seg)];
load(full_model_file);

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


%% Plot raster for all trials from each variant

if ~strcmp(task_variant,'all')
    figure
    for n = 1:n_variants
        subplot(1,n_variants,n)
        trials_plot = task==variants(n);
        modulation_trials = ...
            get_modulation_trials(trial_times(:,trials_plot),...
            infusion_start,infusion_end);
        plot_raster(Y(:,trials_plot), time_ds, [], ...
            modulation_trials(1), modulation_trials(2));
        title([variant_names{n},' trials'])
        yl = ylim;

        if strcmp(task_variant,'rewarded')
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,variants(n));
        else
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,1);
        end
        xlabel('Time since cue onset (seconds)')
        ylabel('Trial #')
    end
    sgtitle('All trials')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot raster for matched trials from each variant
    % Plot the same number of trials for each variant
    raster = plot_raster_variants(Y, time_ds, trial_times, ...
        infusion_start, infusion_end, variants, task);

    for n = 1:n_variants
        subplot(1,n_variants,n)
        hold on
        yl = ylim;
        title([variant_names{n},' trials'])
        if strcmp(task_variant,'rewarded')
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,variants(n));
        else
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,1);
        end
        xlabel('Time since cue onset (seconds)')
        if n == 1
            ylabel('Trial #')
        end
    end
    sgtitle('Equal trials') 
    % Trials from the task variant with fewest trials are matched to 
    % trials from the other variants based on proximity in time.
end

%% Infer spike rates from spike data and estimated model parameters
% Use all trials to infer spike rate
% (note that xK and WK from the estimated model are estimated from 
% half of the trials)

xK_inf = cell(n_models,1);
WK_inf = cell(n_models,1);
Wku_inf = cell(n_models,1);
logll_inf = cell(n_models,1);

for m = 1:n_models
    model_ind = m;
    Y_m = Y(:,trial_mat(:,model_ind))';
    K=size(Y,2);
    
    spikes_task = Y_m;
    dN = zeros(length(complete_trials),size(spikes_all_ms,1));
    dN(trial_mat(:,model_ind),:) = spikes_task;
    missing_trials = ~trial_mat(:,model_ind);
    
    xK0=xK{m}(:,1);

    [xK_inf{m},WK_inf{m},Wku_inf{m},logll_inf{m},SumXkTerms,sumPPll]= ...
        DecodingAlgorithms_IG.PPSS_EStep_missing(A,Qhat{m},xK0,dN,...
        missing_trials,Hk,fitType,delta,gammahat{m},n_seg);
end

%% Plot estimated spike rates (not including history effect)
maxrate = 0;
minrate = 0;
for n = 1:n_variants
    maxrate = max([maxrate,max(max(exp(xK_inf{n})/delta))]);
    minrate = min([minrate,min(min(exp(xK_inf{n})/delta))]);
end

t_ds = (1:n_seg)*l_seg-l_seg/2;

figure
for n = 1:n_variants
    subplot(1,n_variants,n)
    imagesc(t_ds,mean(trial_times)/60,...
        [exp(xK_inf{n})/delta]');
    c=colorbar;
    axis xy;
    %caxis([0,40])
    caxis([0.1,maxrate])
    xlabel('Time (s)')
    ylabel('Trial time (min)')
    ylabel(c,'Firing rate (spikes/second)')
    cmap = load('glm_cmap');
    cmap = colormap(gca,cmap.cmap);
    set(gca,'ColorScale','log')
    title([variant_names{n},' trials'])
end


%% Compute and plot task encoding effect 
% Plot periods: baseline, inhibition, and recovery

before_trials = (trial_times(2,:)<infusion_start(1))';
after_trials = (trial_times(1,:)>infusion_end(end))';

n_mc = 1000;

last_before = find(before_trials,1,'last');
first_after = find(after_trials,1,'first');
last_after = n_trials;

trials_per_seg = 20;
trial_segs = zeros(3,trials_per_seg);

n_ds = 50;
srf_ci = cell(n_models,3);

figure
ax = [];
for m = 1:n_models
    subplot(1,n_models,m)
    model_trials = find(trial_mat(:,model_ind));
    trials_seg(1,:) = ...
        find(model_trials<=last_before,trials_per_seg,'last');
    trials_seg(2,:) = ...
        find(model_trials>=first_after,trials_per_seg,'first');
    trials_seg(3,:) = ...
        find(model_trials<=last_after,trials_per_seg,'last');
    recovery_time = mean(trial_times(:,model_trials(trials_seg(3,:))));
    for l = 1:3
        trials_run = model_trials(trials_seg(l,:));
        srf = zeros(n_s/n_ds,trials_per_seg);
        srf_mc_j = zeros(n_s/n_ds,n_mc*trials_per_seg);
        srf_mc = zeros(n_s/n_ds,n_mc);
        srf_ci{m,l} = zeros(n_s/n_ds,2);
        for j = 1:trials_per_seg
            trial_run = trials_run(j);
            k=trial_run;
            trial_ind = trial_run;


            Yk = Y(:,trial_ind)';
            nst = nspikeTrain( (find(Yk==1)-1)*delta);
            nst.setMinTime(minTime);
            nst.setMaxTime(maxTime);
            Hk_trial = histObj.computeHistory(nst).dataToMatrix;
            
            histEffect=exp(gammahat{m}*Hk_trial')';
            stimK=basisMat*xK_inf{m}(:,k);
            stimEffect=exp(stimK);
            mu = histEffect.*stimEffect;
            stimEffect_scaled = stimEffect/mean(mu);
            mu = histEffect.*stimEffect_scaled;
            mu = reshape(mu,n_ds,[]);
            srf(:,j) = (mean(mu)); 
            
            % Estimate confidence intervals with Monte Carlo
            for n = 1:n_mc
                R = mvnrnd(xK_inf{m}(:,k),WK_inf{m}(:,:,k),1);
                stimK=basisMat*R';
                stimEffect=exp(stimK);
                mu = histEffect.*stimEffect;
                stimEffect_scaled = stimEffect/mean(mu);
                mu = histEffect.*stimEffect_scaled;
                mu = reshape(mu,n_ds,[]);
                srf_mc_j(:,n+(j-1)*n_mc) = (mean(mu));
            end
        end
        for n = 1:n_mc
            ind = [];
            for j = 1:length(trials_run)
                ind = [ind, n+(j-1)*n_mc];
            end
            srf_mc(:,n) = mean(srf_mc_j(:,ind),2);
        end

        srf = mean(srf,2)';
        for n = 1:(n_s/n_ds)
            srf_ci{m,l}(n,1) = quantile(srf_mc(n,:),0.025);
            srf_ci{m,l}(n,2) = quantile(srf_mc(n,:),0.975);    
        end
        time_plot = mean(reshape(time_ds,n_ds,[]));
        plotConfInterval_cont(time_plot,...
            srf,srf_ci{m,l}(:,1)', srf_ci{m,l}(:,2)',l-1)
        hold on
    end
    ax = [ax,gca];

    xlabel('Time (s)')
    ylabel('Firing rate modulation');
    title([variant_names{m},' trials'])
end

linkaxes(ax,'xy')
for m = 1:n_models
    subplot(1,n_models,m)
    yl = ylim;
    plot_trial_info(yl,sample_start,sample_end,...
        reward_start,reward_end,match_start,0.2,m);
end