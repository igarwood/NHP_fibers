function lfp_data = LFP_analysis(session, suppress_figs, remove_noise)

%% Parameters
% Parameters for notch filter:
filter_order = 2;
filter_lowpass = 59.6; % Hz
filter_highpass = 60.4; % Hz

% Parameters for multitaper analysis
moving_win = [30,10];
tapers = [10,19];

% Upper LFP
upper_LFP = 100;

spec_freq_max = 100; % Hz, for plotting
plot_unit = 'min'; % 'min' or 'sec'

% If there are periods of high amplitude noise, choose whether to remove
% data including and following noise with remove_noise parameter
noise_threshold = 10; % defined empirically

%% Extract session information:
[session_info] = setup(session);

folder = session_info.folder;
fileroot = session_info.fileroot;
pharmacology = session_info.pharmacology;
fs = session_info.fs_LFP;

file = [folder,fileroot];

load(file,'lfp');
lfp_elec = lfp.lfp_elec;
lfp = lfp.lfp_data;
n_elecs = size(lfp,1);
N = size(lfp,2);
time = (0:(N-1))/fs;


infusion_start = pharmacology.infusion_start;
infusion_end = pharmacology.infusion_end;
infusion_rate = pharmacology.infusion_rate;

%% Filter lfp

d = designfilt('bandstopiir','FilterOrder',filter_order, ...
               'HalfPowerFrequency1',filter_lowpass,...
               'HalfPowerFrequency2',filter_highpass, ...
               'DesignMethod','butter','SampleRate',fs);

for i = 1:n_elecs
    lfp(i,:) = filtfilt(d,lfp(i,:));
end 

%% Calculate spectrogram:
[time,spec_time, spec_freq, lfp_chan, spectrogram] = ...
    spectData(lfp(lfp_elec,:),'time',time,'fs',fs,...
    'moving_win',moving_win,'tapers',tapers);

spectrogram = spectrogram/((upper_LFP/1.96)^2/fs);
[~,spec_freq_max_ind] = min(abs(spec_freq-spec_freq_max));


% Optionally remove consecutive periods of noise from analysis
if remove_noise
    rms_spectrogram = rms(movmean(spectrogram(:,1:spec_freq_max_ind)',...
        10,2));
    noise = rms_spectrogram>noise_threshold;
    spectrogram = spectrogram(1:find(noise,1),:);
    spec_time = spec_time(1:find(noise,1));
end

%% Plots
if ~suppress_figs
    if strcmp(plot_unit,'min')
        time_scale = 60;
    else
        time_scale = 1;
    end

    [~,spec_freq_max_ind] = min(abs(spec_freq-spec_freq_max));
    figure

    plotSpectrogram(spectrogram(:,1:spec_freq_max_ind), ...
        spec_time/time_scale,spec_freq(1:spec_freq_max_ind));
    hold on 
    yl = [0,spec_freq_max];
    ylim(yl);

    for n = 1:length(infusion_start)
        inject_times = [(infusion_start(n));infusion_end(n)]/time_scale;
        plot(inject_times,ones(size(inject_times))*(yl(2)+1),'k',...
            'linewidth',5)
    end

    ylim([0,yl(2)+2])
    %caxis([-10,40])
    xlabel(['Time (',plot_unit,')'])
    ax = gca;

    figure
    plot(time/time_scale,lfp(lfp_elec,:),'k')
    ylabel('Amplitude (mV)')
    box off
    ax = [ax,gca];


    linkaxes(ax,'x');
    xlim([0,max(spec_time/time_scale)])
end

%% Output structure
lfp_data.lfp = lfp;
lfp_data.spectrogram = spectrogram;
lfp_data.time = time;
lfp_data.spec_time = spec_time;
lfp_data.spec_freq = spec_freq;
if ~suppress_figs
    lfp_data.ax = ax;
    lfp_data.plot_unit = plot_unit;
end

end