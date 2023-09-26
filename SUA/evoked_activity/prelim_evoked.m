%% Preliminary analysis of evoked activity 
% Plot rasters and compute psth before and after intracranial injections


session = 'pmc_gaba2';
spike_ind = 1;
task_variant = 'rewarded';
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
end

%% Extract session information:
[session_info] = setup(session);
folder = session_info.folder;
fileroot = session_info.fileroot;
pharmacology = session_info.pharmacology;

infusion_start = pharmacology.infusion_start;
infusion_end = pharmacology.infusion_end;
infusion_rate = pharmacology.infusion_rate;

fs = session_info.fs_spike;

load('trial_info.mat')

%% Extract task information:
[spikes_all, time, task, trials, trial_times, correct_trials, ...
    complete_trials] = evoked_spike(session,spike_ind,task_variant);

trial_dur = time(trials(2,1))-time(trials(1,1));
trial_time = time(trials(1,1):trials(2,1))-time(trials(1,1));

%% Remove outlier trials
% Noise can be improperly sorted into units. Remove trials corrupted by
% noise

for n = 1:2
    trial_outliers = remove_spike_outliers(spikes_all,0.12*fs);
    spikes_all(:,trial_outliers) = [];
    task(trial_outliers) = [];
    trials(:,trial_outliers) = [];
    trial_times(:,trial_outliers) = [];
    complete_trials(trial_outliers) = [];
    correct_trials(trial_outliers) = [];
    n_trials = size(spikes_all,2);
end

%% Plot rasters
% Plots 3 types of rasters (see function), [2+n_variants] figures total
y_axis_time = 0; % change to 1 for y_axis in minutes instead of trial #
raster = plot_all_rasters(spikes_all,trial_time, trial_times, ...
    infusion_start, infusion_end, task_variant, variants, variant_names,...
    task, y_axis_time);

%% Plot PSTH before/after infusions
psth_win = 0.1*fs; % number of samples per window
psth_res = 0.05*fs; % step size between windows

before_trials = trial_times(1,:) < infusion_start(1);
after_trials = trial_times(2,:) > infusion_end(end);

spikes_before = spikes_all(:,before_trials);
spikes_after = spikes_all(:,after_trials);
task_before = task(before_trials);
task_after = task(after_trials);

before_psth = zeros(n_variants,length(trial_time)-1);
after_psth = zeros(n_variants,length(trial_time)-1);
figure
for n = 1:n_variants
    before_spikes = spikes_before(:,task_before==variants(n));
    n_trial_var = size(before_spikes,2);
    before_spikes_avg = sum(before_spikes,2)/n_trial_var;
    before_psth(n,:) = movmean(before_spikes_avg,psth_win)*fs; 
    
    after_spikes = spikes_after(:,task_after==variants(n));
    n_trial_var = size(after_spikes,2);
    after_spikes_avg = sum(after_spikes,2)/n_trial_var;
    after_psth(n,:) = movmean(after_spikes_avg,psth_win)*fs; 
    
    subplot(2,n_variants,n)
    plot(trial_time(1:psth_res:end-1),before_psth(n,1:psth_res:end),'k')
    hold on
    title([variant_names{n},' trials'])
    
    
    subplot(2,n_variants, n_variants+n)
    plot(trial_time(1:psth_res:end-1),after_psth(n,1:psth_res:end),'k')
    hold on
    xlabel('Time since cue onset')
    
    for m = 1:2
        subplot(2,n_variants,n_variants*(m-1)+n)
        yl = ylim;
        if strcmp(task_variant,'rewarded')
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,variants(n));
        else
            plot_trial_info(yl,sample_start,sample_end,...
                reward_start,reward_end,match_start,0.2,1);
        end
        if n == 1
            ylabel('Firing rate (spikes/second)')
        end
    end
    
end

