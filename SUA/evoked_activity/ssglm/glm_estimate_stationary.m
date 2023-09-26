%% Estimate stationary generalized linear model from baseline trials
% Note that this uses nSTAT toolbox

session = 'pmc_gaba2';
spike_ind = 1;
task_variant = 'rewarded';
n_seg = 50;
trial_dur = 6;

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
before_trials = (trial_times(2,:)<infusion_start(1))';
trials_keep = before_trials; % Change if any trials are not to be included
test_trials = 2:2:sum(before_trials); % Trials to be held out from estimation

n_trials = sum(trials_keep);
Y = Y(:,trials_keep);
trial_times = trial_times(:,trials_keep);
complete_trials = complete_trials(trials_keep);
task = task(trials_keep);
complete_trials(test_trials) = 0; 
% ^ ss-glm will be estimated from complete trials

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


%% Plot raster for all trials
figure
plot_raster(Y, time_ds);
title('All trials')
yl = ylim;
plot_trial_info(yl,sample_start,sample_end,...
    reward_start,reward_end,match_start,0.2,1)
xlabel('Time since cue onset (seconds)')
ylabel('Trial #')

%% Plot raster for all trials from each variant
if ~strcmp(task_variant,'all')
    n_variants = length(variants);

    for n = 1:n_variants
        figure
        trials_plot = task==variants(n);
        plot_raster(Y(:,trials_plot), time_ds);
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

end
%% GLM data structures:
xK = cell(n_models,1);
WK = cell(n_models,1);
WkuFinal = cell(n_models,1);
Qhat = cell(n_models,1);
gammahat = cell(n_models,1);
fitResults = cell(n_models,1);
stimulus = cell(n_models,1);
stimCIs = cell(n_models,1);
logll = cell(n_models,1);
QhatAll = cell(n_models,1);
gammahatAll = cell(n_models,1);
nIter = cell(n_models,1);
Q0 = cell(n_models,1);
fR = cell(n_models,1);


%% Estimate GLM
dev = zeros(n_models,1);

est_var_yn = 0;
for m = 1:n_models
    model_ind = m;
    n_trials_m = sum(trial_mat(:,model_ind));
    B_init = zeros(R+J,n_models);
    
    init_inds = find(trial_mat(:,model_ind));
    
    Y_flat = reshape(Y(:,init_inds),[],1);
    spike_hist_flat = reshape(spike_hist(:,init_inds,:),[],J);
    trial_segs_flat = repmat(basisMat,length(init_inds),1);

    X_flat = [trial_segs_flat,spike_hist_flat];
    tic
    [B_init(:,model_ind),dev(model_ind),...
        stats{model_ind}] = glmfit(X_flat,Y_flat,...
        'poisson','constant','off','link','log');
    t=toc;
    disp(['Baseline model ',num2str(m),' estimated in ',num2str(t), ...
        ' seconds'])
    
    xK{m} = B_init(1:R,model_ind);
    gammahat{m} = B_init(end-J+1:end,model_ind)';


end
   
save([save_folder,session,'_unit',num2str(spike_ind),'_',task_variant,...
    '_','R',num2str(R),'_stationary'],'xK','gammahat','stats','dev');


%% Plot estimated spike rates (not including history effect)
maxrate = 0;
minrate = 0;
for n = 1:n_variants
    maxrate = max([maxrate,max(max(exp(xK{n})/delta))]);
    minrate = min([minrate,min(min(exp(xK{n})/delta))]);
end

t_ds = (1:n_seg)*l_seg-l_seg/2;


figure
for n = 1:n_variants
    subplot(1,n_variants,n)
    imagesc(t_ds,mean(trial_times)/60,...
        [exp(xK{n})/delta]');
    title([variant_names{n},' trials'])
    c=colorbar;
    axis xy;
    caxis([0.1,maxrate])
    xlabel('Time (s)')
    ylabel('Trial time (min)')
    ylabel(c,'Firing rate (spikes/second)')
    cmap = load('glm_cmap');
    cmap = colormap(gca,cmap.cmap);
    %set(gca,'ColorScale','log')
end

