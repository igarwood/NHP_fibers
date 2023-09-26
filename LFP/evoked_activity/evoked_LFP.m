function [lfp_all, time, task, trials, trial_times, correct_trials,...
    complete_trials] = evoked_LFP(session, task_variant)
% Options for task_variant:
% 'rewarded' for rewarded v non-rewarded (aka correct v incorrect) trials
% 'choice' for the choice made during the match phase (1, 2, or 3);
% 'sample' for the sample shown during the sample phase
%       - only include correct trials
% 'lr' for left v right saccade location 
% If no task_variant input, task is vector of zeros
if nargin == 1
    task_variant = [];
end

%% Extract session information:
[session_info] = setup(session);
folder = session_info.folder;
fileroot = session_info.fileroot;
pharmacology = session_info.pharmacology;
fs = session_info.fs_spike;

file = [folder,fileroot];

infusion_start = pharmacology.infusion_start;
infusion_end = pharmacology.infusion_end;
infusion_rate = pharmacology.infusion_rate;

%% Load lfp and task information:
load(file,'lfp','fs_lfp','task_info')

lfp = lfp.lfp_data(lfp.lfp_elec,:);
time = (0:(length(lfp)-1))/fs_lfp;


%% Extract task information:
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


trial_dur = fs_lfp*6;

%% Determine trial times and indices
trial_times = (trials-1)/fs_task;
n_trials = size(trials,2);

%% Assign LFP to trials
lfp_all = zeros(trial_dur,n_trials);
ind_offset = 0;

for n = 1:n_trials
    trial_time = trial_times(:,n);
    trial_locs = round(trial_time*fs_lfp)+1;

    lfp_all(:,n) = lfp(trial_locs(1):(trial_locs(2)-1));
end

end