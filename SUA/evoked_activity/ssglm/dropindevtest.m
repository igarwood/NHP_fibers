%% Drop-in-deviance test aka LR test. 
%
% Dependencies: NHP_fibers and nSTAT toolboxes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To determine if there is significant task encoding, first estimate a full 
% and reduced model from all trials (using ssglm_estimate)
% Add those files + corresponding information here:

task_variant = 'all';

spike_ind = 1;
session = 'pmc_gaba2';
numBasis_full = 50;
numBasis_red = 1;
n_history = 9;
trial_dur = 6;


[session_info] = setup(session);
save_folder = session_info.save_folder;
full_model_file = [save_folder,session,'_unit',num2str(spike_ind),'_',...
    task_variant,'_','R',num2str(numBasis_full)];
reduced_model_file = [save_folder,session,'_unit',num2str(spike_ind),...
    '_',task_variant,'_','R',num2str(numBasis_red)];

%% Load spike data
[spikes_all, time, task, trials, trial_times, correct_trials, ...
    complete_trials] = evoked_spike(session,spike_ind,task_variant);

%% Extract trial data and compute spike history:

delta = 1e-3;
sampleRate=1/delta;
fitType = 'poisson';


n_trials = length(correct_trials);
spikes_all_ms = zeros(trial_dur*1000,n_trials);

fs = 30000;
s_time = 0:1/fs:(length(spikes_all)-1)/fs;
t_ds = 0.001;
s_ds = fs*t_ds;
n_s = trial_dur/t_ds;

time_ds = zeros(n_s,1);

for n = 1:n_s
    ind = ((n-1)*s_ds+1):(n*s_ds);
    time_ds(n) = mean(s_time(ind));
    spikes_all_ms(n,:) = sum(spikes_all(ind,:));
    spikes_all_ms(n,spikes_all_ms(n,:)>1) = 1;
end
Y = spikes_all_ms;
for n = 1:2
    trial_outliers = remove_spike_outliers(Y,120);
    Y(:,trial_outliers) = [];
    trial_times(:,trial_outliers) = [];
    complete_trials(trial_outliers) = [];
    correct_trials(trial_outliers) = [];
    n_trials = size(Y,2);
end

trials_keep = 1:n_trials;
trials_remove = 1:2:n_trials;

n_trials = length(trials_keep);
Y = Y(:,trials_keep);
trial_times = trial_times(:,trials_keep);
complete_trials = ones(size(trials_keep));
complete_trials(trials_remove) = 0;
model_trials = complete_trials;
missing_trials = ~model_trials;
Y_m = Y(:,complete_trials==1)';
dN = zeros(length(complete_trials),size(Y,1));
dN(complete_trials==1,:) = Y_m;
windowTimes = [0,1,5,10,15,20,30,40,50,100]/1000;

K = size(Y,2);
Hk = cell(1,K);
Hk_null = cell(1,K);
minTime = 0;
maxTime=(size(Y_m,2)-1)*delta;
histObj = History(windowTimes,minTime,maxTime);
for k=1:K
    nst{k} = nspikeTrain( (find(Y(:,k)==1)-1)*delta);
    nst{k}.setMinTime(minTime);
    nst{k}.setMaxTime(maxTime);
    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
end


%% Drop in dev test:

basisWidth = (maxTime-minTime)/numBasis_full;
unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,...
    minTime,maxTime,sampleRate);
basisMat_full = unitPulseBasis.data;

basisWidth = (maxTime-minTime)/numBasis_red;
unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,...
    minTime,maxTime,sampleRate);
basisMat_red = unitPulseBasis.data;

full_model = load(full_model_file);
red_model = load(reduced_model_file);        
        
A=eye(numBasis_full,numBasis_full);
x0=full_model.xK{1}(:,1);
Q = full_model.Qhat{1};
gamma = full_model.gammahat{1};

[xK,~,~,logll_full,~,~]...
    =DecodingAlgorithms_IG.PPSS_EStep_missing(A,Q,x0,dN,...
    missing_trials,Hk,fitType,delta,gamma,numBasis_full);

A=eye(numBasis_red,numBasis_red);
x0=red_model.xK{1}(:,1);
Q = red_model.Qhat{1};
gamma = red_model.gammahat{1};
[~,~,~,logll_red,~,~]...
    =DecodingAlgorithms_IG.PPSS_EStep_missing(A,Q,x0,dN,...
    missing_trials,Hk,fitType,delta,gamma,numBasis_red);


% number of params = (number of trials + 1(Q0))*num_basis + numhistory
n_trials = size(full_model.xK{1},2);
param_full = (n_trials+1)*numBasis_full+n_history;
param_red = (n_trials+1)*numBasis_red+n_history;
param_diff = param_full-param_red;
dev_full = -2*logll_full;
dev_red = -2*logll_red;

drop_in_dev = dev_red-dev_full;

p_val = 1-chi2cdf(drop_in_dev,param_diff);

fprintf(['p-value for drop in deviance test = ', num2str(p_val) '\n'])



