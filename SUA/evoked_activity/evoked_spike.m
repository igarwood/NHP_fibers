function [spikes_all, time, task, trials, trial_times, correct_trials,...
    complete_trials] = evoked_spike(session, spike_ind, task_variant)
% Options for task_variant:
% 'rewarded' for rewarded v non-rewarded (aka correct v incorrect) trials
% 'choice' for the choice made during the match phase (1, 2, or 3);
% 'sample' for the sample shown during the sample phase
%       - only include correct trials
% 'lr' for left v right saccade location 
% If no task_variant input, task is vector of zeros
if nargin == 2
    task_variant = [];
end
%% Parameters
suppress_spike_figs = 1;

%% Extract session information:
[session_info] = setup(session);
folder = session_info.folder;
save_folder = session_info.save_folder;
fileroot = session_info.fileroot;
pharmacology = session_info.pharmacology;
fs = session_info.fs_spike;

file = [folder,fileroot];
spike_folder = [save_folder,'spike_data/'];
spike_file = [spike_folder,'spikes_cluster_',session];

infusion_start = pharmacology.infusion_start;
infusion_end = pharmacology.infusion_end;
infusion_rate = pharmacology.infusion_rate;

[spike_data] = load_sorted_spikes(session,suppress_spike_figs);
spike_locs = spike_data.spike_locs{spike_ind};
time = spike_data.s_time;

%% Extract task information:
load([folder,session],'task_info','fs_spike');

correct_trials = task_info.rewarded_trials;
complete_trials = task_info.complete_trials;
samples = task_info.samples;
fs_task = task_info.fs_task;

if strcmp(task_variant,'rewarded')
    task = task_info.rewarded_trials;
elseif strcmp(task_variant,'choice')
    task = task_info.choice_id;
elseif strcmp(task_variant,'sample')
    sample_id(correct_trials==0) = -1;
    task = task_info.sample_id;
elseif strcmp(task_variant,'lr')
    saccade_loc = task_info.correct_loc==11;
    saccade_loc((correct_trials==0).*(task_info.correct_loc==12)==1) = 1;
    saccade_loc((correct_trials==0).*(task_info.correct_loc==11)==1) = 0;
    task = saccade_loc;
else
    task = zeros(size(correct_trials));
end
    
% Use only complete trials; start_trial = sample-1s, end_trial = sample+5s
task(complete_trials==0) = [];
correct_trials(complete_trials==0) = [];
samples(:,complete_trials==0) = [];
complete_trials(complete_trials==0) = [];


trials = round([samples(1,:)-fs_task;...
    samples(1,:)+(5*fs_task)]);
trial_dur = fs_spike*6;
%% Insert "trials" between completed trials
% Purpose: to have ~ evenly spaced trials throughout the experiment 
% This is done for the purpose of the ss-glm; neural data from these 
% "trials" are not included for analysis

min_trial_interval = min(trials(1,2:end)-trials(2,1:end-1));
while max(trials(1,2:end)-trials(2,1:end-1))...
        >((min_trial_interval)*2+6*fs_task)
    spaced_trials = find((trials(1,2:end)-trials(2,1:end-1))...
        >((min_trial_interval)*2+6*fs_task));
    for n = 1:length(spaced_trials)
        insert_trial = [trials(2,spaced_trials(n))+min_trial_interval;...
            trials(2,spaced_trials(n))+min_trial_interval+6*fs_task];
        trials_insert = [trials(:,1:spaced_trials(n)),insert_trial,...
            trials(:,spaced_trials(n)+1:end)];
        correct_trials_insert = [correct_trials(1:spaced_trials(n));0;...
            correct_trials(spaced_trials(n)+1:end)];
        task_insert = [task(1:spaced_trials(n));0;...
            task(spaced_trials(n)+1:end)];
        complete_trials_insert = [complete_trials(1:spaced_trials(n));0;...
            complete_trials(spaced_trials(n)+1:end)];   
    end
    trials = trials_insert;
    correct_trials = correct_trials_insert;
    task = task_insert;
    complete_trials = complete_trials_insert;
end

%% Determine trial times and indices
trial_times = [time(trials(1,:));time(trials(2,:))];
trial_ind = mean(trials);
n_trials = size(trials,2);

%% Assign spikes to trials
spikes_all = zeros(trial_dur,n_trials);
ind_offset = 0;

for n = 1:n_trials
    trial_locs = zeros(1,2);
    spikes_all_n = zeros(trial_dur,1);

    trial_time = trial_times(:,n);
    trial_locs = round(trial_time*fs_spike);

    trial_spike_locs = spike_locs((spike_locs>=trial_locs(1))...
        .*(spike_locs<=trial_locs(2))==1);
    trial_spike_locs = trial_spike_locs - trial_locs(1) + 1;
    spikes_all_n(trial_spike_locs) = 1;

    spikes_all(:,n) = spikes_all_n(1:trial_dur);
end

end
