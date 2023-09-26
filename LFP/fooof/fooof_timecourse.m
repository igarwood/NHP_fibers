%% Calculate power in fooof derived frequency bands over time
session = 'put_GABA';
session_type = 'GABA';
session_loc = 'putamen'; % Alternative: 'pmc'
suppress_lfp_figs = 0; % Plot spectrogram + time series
suppress_fooof_figs = 0;

session_info = setup(session);
infusion_start = session_info.pharmacology.infusion_start;
infusion_end = session_info.pharmacology.infusion_end;
infusion_rate = session_info.pharmacology.infusion_rate;


% The next two commands take a moment to run (both load + analyze data)
if strcmp(session_loc,'pmc')
    fooof_data = fooof_analysis_bs(session_type,suppress_fooof_figs); 
    lfp_data = LFP_analysis(session, suppress_lfp_figs,1);
else
    fooof_data = fooof_analysis_put(session_type,suppress_fooof_figs); 
    lfp_data = LFP_analysis(session, suppress_lfp_figs,0);
end



% From FOOOF:
band_freqs = fooof_data.band_freqs;
band_diff = fooof_data.mean_band_diff;
comparison = sign(band_diff);

% From LFP analysis:
spectrogram = lfp_data.spectrogram;
spec_time = lfp_data.spec_time;
spec_freq =lfp_data.spec_freq;


% Remove frequencies around line noise from analysis
ind = ones(size(spec_freq));
ind(((spec_freq>59.5).*(spec_freq <60.5))==1)=0;
ind = ind==1;
spec_freq = spec_freq(ind);
spectrogram = spectrogram(:,ind);

%% Compute power in each band
% Remove frequencies around line noise from analysis
ind = ones(size(spec_freq));
ind(((spec_freq>59.5).*(spec_freq <60.5))==1)=0;
ind = ind==1;
spec_freq = spec_freq(ind);
spectrogram = spectrogram(:,ind);


y = band_power(spectrogram,spec_freq, band_freqs);
H = size(band_freqs,1);

y_baseline = y(spec_time<infusion_start(1),:);
y_cdf = zeros(size(y));

for h = 1:H
    mu = mean(y_baseline(:,h));
    sigma = std(y_baseline(:,h));
    y_cdf(:,h) = cdf('normal',y(:,h),mu,sigma);
end

%% Plots
time_scale = 60;
if suppress_lfp_figs
    ax = [];
    plot_unit = 'min';
else
    ax = lfp_data.ax;
    plot_unit = lfp_data.plot_unit;
end

if strcmp(plot_unit,'min')
    time_scale = 60;
else
    time_scale = 1;
end
    
figure
for h = 1:H
    subplot(H,1,h)
    hold on
    if comparison(h)==1
        threshold_crossings = find(y_cdf(:,h)>0.95);
    else
        threshold_crossings = find(y_cdf(:,h)<0.05); 
    end
    width = (spec_time(2)-spec_time(1))/time_scale;
    for n = 1:length(threshold_crossings)
        rectangle('position',...
            [spec_time(threshold_crossings(n))/time_scale,0,width,1],...
            'facecolor',[0,0,0,0.1],'edgecolor','none')
    end
    
    plot(spec_time/time_scale,y_cdf(:,h),'k')
    ax = [ax,gca];
    hold on
    for n = 1:length(infusion_start)
        plot([infusion_start(n),infusion_end(n)]/time_scale,...
            [1.01,1.01],'k')
    end
    if h == H  
        xlabel(['Time (',plot_unit,')'])
    end
end
linkaxes(ax,'x')
xlim([0,max(spec_time/time_scale)])
ylim([0,1.01])