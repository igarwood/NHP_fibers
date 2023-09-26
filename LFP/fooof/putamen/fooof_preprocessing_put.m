%% fooof_preprocessing
%
% Code for preprocessing LFP data for comparison of LFP oscillatory 
% structure before and after intracranial GABA or SALINE infusions
%
% Written to analyze premotor cortex data from Garwood, et al., 2022
%% Specify sessions:
session = 'put_gaba';

%% Parameters
% Save extracted data for fooof analysis?
save_data = 1;

% Parameters for notch filter:
filter_order = 2;
filter_lowpass = 59.6; % Hz
filter_highpass = 60.4; % Hz

% If there are periods of high amplitude noise, choose whether to remove
remove_noise = 1;
noise_threshold = 60; % defined empirically

% Sample parameters:
sample_period = 10*60; % in seconds, period of time before/after injections
sample_length = 30; % in seconds, duration of each LFP sample

% Bootstrap parameters: 
n_samples = 60; % number of samples per bootstrap samples
n_bs = 100; % number of bootstrap samples

% Parameters for multitaper analysis
tapers = [10,19];
n_freq = 3278; 
% ^ number of multitaper frequency samples (will depend on taper 
% parameters and fs; if not known, run LFP_analysis with the same 
% parameters to see the number of freqs)

% Upper LFP
upper_LFP = 100;
%% Extract session information:
% Assumes that folder, save_folder, and sampling rate is the same for all
% sessions

session_times = zeros(2,2);


session_info = setup(session);
folder = session_info.folder;
save_folder = [session_info.folder,'/data_for_fooof/'];
fs = session_info.fs_LFP;

fileroot = session_info.fileroot;


pharmacology = session_info.pharmacology;
session_times = zeros(2,2);
session_times(1,1) = pharmacology.infusion_start(1)-sample_period;
session_times(1,2) = pharmacology.infusion_start(1);
session_times(2,1) = pharmacology.infusion_end(end);
session_times(2,2) = pharmacology.infusion_end(end)+sample_period;


%% Set up data structures:

samp_length = sample_length*fs;

% Create logical matrix for valid time points to calculate spectra from
% (e.g. all timepoints that have t:t+30 sec of data to calculate from)
lfp_bin = zeros(sample_period*fs,2);

samples_total = n_bs*n_samples*2;
spec_sample_session = zeros(n_samples,n_freq);
spec_sample_all = zeros(samples_total,n_freq);

%% Load LFP and assign to sample periods

d = designfilt('bandstopiir','FilterOrder',filter_order, ...
               'HalfPowerFrequency1',filter_lowpass,...
               'HalfPowerFrequency2',filter_highpass, ...
               'DesignMethod','butter','SampleRate',fs);


load([session_info.folder,fileroot],'lfp');
channel = lfp.lfp_elec;
lfp_i = lfp.lfp_data(lfp.lfp_elec,:);
N = size(lfp.lfp_data,2);
time = (0:(N-1))/fs;

lfp_i = filtfilt(d,lfp_i);

lfp = zeros(sample_period*fs,2);
for j = 1:2
    sample_ind = find((time>session_times(j,1)).* ...
        (time<=(session_times(j,2)))==1);
    ss_length = length(sample_ind);
    % Since some sample periods are less than 10 min long, find valid
    % indices session by session
    ind_start = (ind-1)*ss_length+1;
    ind_end = ind_start+ss_length-1;
    lfp_bin(ind_start:(ind_end-samp_length),j) = 1;
    lfp(ind_start:ind_end,j) = lfp_i(sample_ind)';
end

lfp_ind{1} = find(lfp_bin(:,1));
lfp_ind{2} = find(lfp_bin(:,2));
%% Compute spectral samples

params.Fs=fs;
params.tapers=tapers;

for j = 1:2
    n = 1;
    for n = 1:n_bs
        tic
        noise = [];
        parfor m = 1:n_samples
            if remove_noise
                noise = true;
            end
            sample_ind = 0;
            lfp_sample = zeros(samp_length,1);
            while noise % Collect samples until a non-noise sample is found
                sample_ind = randsample(lfp_ind{j},1);
                lfp_sample = ...
                    lfp(sample_ind:sample_ind+samp_length-1,j);
                noise = rms(lfp_sample)>noise_threshold;
            end
            [spec_sample,f] = mtspectrumc(lfp_sample, params);
            spec_sample = spec_sample/((upper_LFP/1.96)^2/fs);
            spec_sample_session(m,:) = spec_sample(1:n_freq);
        end
        ind = (j-1)*n_bs*n_samples + (n-1)*n_samples +...
            [1:n_samples];
        spec_sample_all(ind,:) = spec_sample_session;
        toc
    end
end

sample_ind = randsample(lfp_ind{j},1);
lfp_sample = lfp(sample_ind:sample_ind+samp_length-1,j);
[~,f] = mtspectrumc_param(lfp_sample, params);
freqs = f(1:n_freq);

%% Visualize output:
cmap = load('white_spec_scale');

figure
subplot(1,3,1)
ind_start = 1;
ind_end = n_bs*n_samples;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_all([ind_start:ind_end],:)'));
axis xy;
caxis([-20,30])
colormap(gca,cmap.cmap)
ylim([0,100])
ylabel('Frequency')
xlabel('Samples')
title('Before GABA')
ax = gca;

subplot(1,3,2)
ind_start = n_bs*n_samples+1;
ind_end = n_bs*n_samples*2;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_all([ind_start:ind_end],:)'));
axis xy;
caxis([-20,30])
colormap(gca,cmap.cmap)
ylim([0,100])
ylabel('Frequency')
xlabel('Samples')
title('After GABA')
ax = [ax,gca];
linkaxes(ax,'xy')

% Rough estimate of power changes:
subplot(1,3,3) 
f_ind_plot = (double(freqs<(filter_lowpass-0.8)) + ...
    double(freqs>(filter_highpass+0.8)))==1;
plot(freqs(f_ind_plot),...
    mean(pow2db(spec_sample_all((1:(n_bs*n_samples)),f_ind_plot))))
hold on
plot(freqs(f_ind_plot),mean(pow2db(spec_sample_all...
    (((n_bs*n_samples+1):(n_bs*n_samples*2)),f_ind_plot))))
ylabel('Mean Power (dB)')
xlabel('Frequency (Hz)')
xlim([0,100])
legend({'Before GABA','After GABA'})

%% Save data

writematrix(spec_sample_all,[save_folder, 'spec_sample_',session,'.csv'])
writematrix(freqs',[save_folder,'put_freqs.csv'])
