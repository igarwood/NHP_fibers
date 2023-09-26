%% fooof_preprocessing
%
% Code for preprocessing LFP data for comparison of LFP oscillatory 
% structure before and after intracranial GABA or SALINE infusions
%
% Written to analyze premotor cortex data from Garwood, et al., 2022
%% Specify sessions:
GABA_sessions = {'pmc_gaba1','pmc_gaba2','pmc_gaba3','pmc_gaba4'};
saline_sessions = {'pmc_sal1','pmc_sal2','pmc_sal3'};

%% Parameters
% Save extracted data for fooof analysis?
save_data = 1;

% Parameters for notch filter:
filter_order = 2;
filter_lowpass = 59.6; % Hz
filter_highpass = 60.4; % Hz

% If there are periods of high amplitude noise, choose whether to remove
remove_noise = 1;
noise_threshold = 55; % defined empirically

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

n_GABA = length(GABA_sessions);
n_saline = length(saline_sessions);

GABA_session_times = cell(1,n_GABA);
saline_session_times = cell(1,n_saline);

GABA_fileroot = cell(1,n_GABA);
saline_fileroot = cell(1,n_saline);

GABA_channels = zeros(1,n_GABA);
saline_channels = zeros(1,n_saline);

for n = 1:n_GABA
    session_info = setup([GABA_sessions{n}]);

    folder = session_info.folder;
    save_folder = [session_info.folder,'/data_for_fooof/'];
    fs = session_info.fs_LFP;
    
    GABA_fileroot{n} = session_info.fileroot;
   
    pharmacology = session_info.pharmacology;
    session_times = zeros(2,2);
    session_times(1,1) = pharmacology.infusion_start(1)-sample_period;
    session_times(1,2) = pharmacology.infusion_start(1);
    session_times(2,1) = pharmacology.infusion_end(end);
    session_times(2,2) = pharmacology.infusion_end(end)+sample_period;
    
    GABA_session_times{n} = session_times;
end

for n = 1:n_saline
    [session_info] = setup([saline_sessions{n}]);
    saline_fileroot{n} = session_info.fileroot;
    
    pharmacology = session_info.pharmacology;
    session_times = zeros(2,2);
    session_times(1,1) = pharmacology.infusion_start(1)-sample_period;
    session_times(1,2) = pharmacology.infusion_start(1);
    session_times(2,1) = pharmacology.infusion_end(end);
    session_times(2,2) = pharmacology.infusion_end(end)+sample_period;
    
    saline_session_times{n} = session_times;
end


%% Set up data structures:

samp_length = sample_length*fs;
lfp_GABA = zeros(sample_period*fs*n_GABA,2);
% Create logical matrix for valid time points to calculate spectra from
% (e.g. all timepoints that have t:t+30 sec of data to calculate from)
lfp_GABA_bin = zeros(sample_period*fs*n_GABA,2);
lfp_GABA_session_tags = zeros(sample_period*fs*n_GABA,1);
for n = 1:n_GABA
    lfp_GABA_session_tags(((n-1)*sample_period*fs+1):...
        (n*sample_period*fs)) = ones(sample_period*fs,1)*n;
end

lfp_saline = zeros(sample_period*fs,2);
% Create logical matrix for valid time points to calculate spectra from
% (e.g. all timepoints that have t:t+30 sec of data to calculate from)
lfp_saline_bin = zeros(sample_period*fs*n_saline,2);
lfp_saline_session_tags = zeros(sample_period*fs*n_saline,1);
for n = 1:n_saline
    lfp_saline_session_tags(((n-1)*sample_period*fs+1):...
        (n*sample_period*fs)) = ones(sample_period*fs,1)*n;
end

samples_total_GABA = n_bs*n_samples*n_GABA*2;
samples_total_saline = n_bs*n_samples*n_saline*2;

spec_sample_GABA = zeros(samples_total_GABA,n_freq);
spec_sample_saline = zeros(samples_total_saline,n_freq);

session_tags_GABA = zeros(samples_total_GABA,1);
session_tags_saline = zeros(samples_total_saline,1);

spec_sample_GABA_session = zeros(n_samples*n_GABA,n_freq);
session_tags_GABA_session = zeros(n_samples*n_GABA,1);

spec_sample_saline_session = zeros(n_samples*n_saline,n_freq);
session_tags_saline_session = zeros(n_samples*n_saline,1);

%% Load LFP and assign to sample periods

d = designfilt('bandstopiir','FilterOrder',filter_order, ...
               'HalfPowerFrequency1',filter_lowpass,...
               'HalfPowerFrequency2',filter_highpass, ...
               'DesignMethod','butter','SampleRate',fs);


for ind = 1:n_GABA
    fileroot = GABA_fileroot{ind};
    
    load([session_info.folder,fileroot],'lfp');
    GABA_channels(ind) = lfp.lfp_elec;
    lfp_i = lfp.lfp_data(lfp.lfp_elec,:);
    N = size(lfp.lfp_data,2);
    time = (0:(N-1))/fs;
    
    channel = GABA_channels(ind);

    lfp_i = filtfilt(d,lfp_i);
    for j = 1:2
        sample_ind = find((time>GABA_session_times{ind}(j,1)).* ...
            (time<=(GABA_session_times{ind}(j,2)))==1);
        ss_length = length(sample_ind);
        % Since some sample periods are less than 10 min long, find valid
        % indices session by session
        ind_start = (ind-1)*ss_length+1;
        ind_end = ind_start+ss_length-1;
        lfp_GABA_bin(ind_start:(ind_end-samp_length),j) = 1;
        lfp_GABA(ind_start:ind_end,j) = lfp_i(sample_ind)';
    end
end

for ind = 1:n_saline
    fileroot = saline_fileroot{ind};
    
    load([session_info.folder,fileroot],'lfp');
    saline_channels(ind) = lfp.lfp_elec;
    lfp_i = lfp.lfp_data(lfp.lfp_elec,:);
    N = size(lfp,2);
    time = (0:(N-1))/fs;
    
    channel = saline_channels(ind);
    
    lfp_i = filtfilt(d,lfp_i);
    for j = 1:2
        sample_ind = find((time>saline_session_times{ind}(j,1)).* ...
            (time<=(saline_session_times{ind}(j,2)))==1);
        ss_length = length(sample_ind);
        % Since some sample periods are less than 10 min long, find valid
        % indices session by session
        ind_start = (ind-1)*ss_length+1;
        ind_end = ind_start+ss_length-1;
        lfp_saline_bin(ind_start:(ind_end-samp_length),j) = 1;
        lfp_saline(ind_start:ind_end,j) = lfp_i(sample_ind)';        
    end
end

lfp_GABA_ind{1} = find(lfp_GABA_bin(:,1));
lfp_GABA_ind{2} = find(lfp_GABA_bin(:,2));
lfp_saline_ind{1} = find(lfp_saline_bin(:,1));
lfp_saline_ind{2} = find(lfp_saline_bin(:,2));
%% Compute spectral samples

params.Fs=fs;
params.tapers=tapers;

for j = 1:2
    n = 1;
    for n = 1:n_bs
        tic
        noise = [];
        parfor m = 1:(n_samples*n_GABA)
            if remove_noise
                noise = true;
            end
            sample_ind = 0;
            lfp_sample = zeros(samp_length,1);
            while noise % Collect samples until a non-noise sample is found
                sample_ind = randsample(lfp_GABA_ind{j},1);
                lfp_sample = ...
                    lfp_GABA(sample_ind:sample_ind+samp_length-1,j);
                noise = rms(lfp_sample)>noise_threshold;
            end
            [spec_sample,f] = mtspectrumc(lfp_sample, params);
            spec_sample = spec_sample/((upper_LFP/1.96)^2/fs);
            spec_sample_GABA_session(m,:) = spec_sample(1:n_freq);
            session_tags_GABA_session(m) = ...
                lfp_GABA_session_tags(sample_ind);
        end
        ind = (j-1)*n_bs*n_samples*n_GABA + (n-1)*n_samples*n_GABA +...
            [1:(n_samples*n_GABA)];
        spec_sample_GABA(ind,:) = spec_sample_GABA_session;
        session_tags_GABA(ind) = session_tags_GABA_session;
        toc
    end
end

for j = 1:2
    n = 1;
    for n = 1:n_bs
        tic
        noise = [];
        parfor m = 1:(n_samples*n_saline)
            if remove_noise
                noise = true;
            end
            sample_ind = 0;
            lfp_sample = zeros(samp_length,1);
            while noise % Collect samples until a non-noise sample is found
                sample_ind = randsample(lfp_saline_ind{j},1);
                lfp_sample = ...
                    lfp_saline(sample_ind:sample_ind+samp_length-1,j);
                noise = rms(lfp_sample)>noise_threshold;
            end
            [spec_sample,f] = mtspectrumc(lfp_sample, params);
            spec_sample = spec_sample/((upper_LFP/1.96)^2/fs);
            spec_sample_saline_session(m,:) = spec_sample(1:n_freq);
            session_tags_saline_session(m) = ...
                lfp_saline_session_tags(sample_ind);
        end
        ind = (j-1)*n_bs*n_samples*n_saline + (n-1)*n_samples*n_saline +...
            [1:(n_samples*n_saline)];
        spec_sample_saline(ind,:) = spec_sample_saline_session;
        session_tags_saline(ind) = session_tags_saline_session;
        toc
    end
end

sample_ind = randsample(lfp_GABA_ind{j},1);
lfp_sample = lfp_GABA(sample_ind:sample_ind+samp_length-1,j);
[~,f] = mtspectrumc_param(lfp_sample, params);
freqs = f(1:n_freq);

%% Visualize output:
cmap = load('white_spec_scale');

figure
subplot(1,4,1)
ind_start = 1;
ind_end = n_bs*n_samples*n_GABA;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_GABA([ind_start:ind_end],:)'));
axis xy;
caxis([-20,30])
colormap(gca,cmap.cmap)
ylim([0,100])
ylabel('Frequency')
xlabel('Samples')
title('Before GABA')
ax = gca;

subplot(1,4,2)
ind_start = n_bs*n_samples*n_GABA+1;
ind_end = n_bs*n_samples*n_GABA*2;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_GABA([ind_start:ind_end],:)'));
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
subplot(1,4,3) 
f_ind_plot = (double(freqs<(filter_lowpass-0.8)) + ...
    double(freqs>(filter_highpass+0.8)))==1;
plot(freqs(f_ind_plot),...
    mean(pow2db(spec_sample_GABA((1:(n_bs*n_samples*n_GABA)),f_ind_plot))))
hold on
plot(freqs(f_ind_plot),mean(pow2db(spec_sample_GABA...
    (((n_bs*n_samples*n_GABA+1):(n_bs*n_samples*n_GABA*2)),f_ind_plot))))
ylabel('Mean Power (dB)')
xlabel('Frequency (Hz)')
xlim([0,100])
legend({'Before GABA','After GABA'})

% Verify that each session is ~equally represented:
subplot(1,4,4)
histogram(session_tags_GABA,'normalization','probability');
xlabel('Session Number')
ylabel('Proportion of samples')

figure
subplot(1,4,1)
ind_start = 1;
ind_end = n_bs*n_samples*n_saline;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_saline([ind_start:ind_end],:)'));
axis xy;
caxis([-20,30])
colormap(gca,cmap.cmap)
ylim([0,100])
ylabel('Frequency')
xlabel('Samples')
title('Before saline')
ax = gca;

subplot(1,4,2)
ind_start = n_bs*n_samples*n_saline+1;
ind_end = n_bs*n_samples*n_saline*2;
imagesc([1:(ind_end-ind_start+1)],freqs,...
    pow2db(spec_sample_saline([ind_start:ind_end],:)'));
axis xy;
caxis([-20,30])
colormap(gca,cmap.cmap)
ylim([0,100])
ylabel('Frequency')
xlabel('Samples')
title('After saline')
ax = [ax,gca];
linkaxes(ax,'xy')

% Rough estimate of power changes:
subplot(1,4,3)
f_ind_plot = (double(freqs<(filter_lowpass-0.8)) + ...
    double(freqs>(filter_highpass+0.8)))==1;
plot(freqs(f_ind_plot),...
    mean(pow2db(spec_sample_saline((1:(n_bs*n_samples*n_saline)),...
    f_ind_plot))))
hold on
plot(freqs(f_ind_plot),mean(pow2db(spec_sample_saline...
    (((n_bs*n_samples*n_saline+1):(n_bs*n_samples*n_saline*2)),...
    f_ind_plot))))
ylabel('Mean Power (dB)')
xlabel('Frequency (Hz)')
xlim([0,100])
legend({'Before saline','After saline'})

% Verify that each session is ~equally represented:
subplot(1,4,4)
histogram(session_tags_saline,'normalization','probability');
xlabel('Session Number')
ylabel('Proportion of samples')
%% Save data

writematrix(spec_sample_GABA,[save_folder, 'spec_sample_GABA.csv'])
writematrix(spec_sample_saline,[save_folder, 'spec_sample_saline.csv'])
writematrix(freqs',[save_folder,'freqs.csv'])
save([save_folder, 'session_tags'],'session_tags_saline',...
    'session_tags_GABA')