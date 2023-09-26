function fooof_data = fooof_analysis_put(sessions,suppress_figs)
%% Parameters
% Parameters from FOOOF preprocessing:
n_samples = 60; % number of samples per bootstrap samples
n_bs = 100; % number of bootstrap samples

% Parameters from FOOOF analysis:
max_peaks = 6; 

% Post-processing parameters:
n_bands = 7; % Number of frequency bands identified; empirically set 
n_runs_km = 10; % Number of runs for kmeans clustering of freq bands
conf_level = 95; % Confidence level to use for confidence intervals

% Default colors for plots:
defaults = [0.6350    0.0780    0.1840;...
        0.0000    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330];
    
%% Get file info:
info = setup;
data_folder = [info.folder,'/fdata_for_fooof/'];
fileroot = 'fooof_output';

%% Load data

GABA = load([folder,fileroot,'_put_GABA_',num2str(max_peaks),...
    '.mat']);
noic = load([folder,fileroot,'_put_noic_',num2str(max_peaks),...
    '.mat']);

spec_freq = readmatrix([folder,'freqs.csv']);
load([folder,'session_tags_put_GABA']);

n_GABA = size(GABA.peak_cf,1);

%% Determine frequency bands of interest:

all_cf_GABA = reshape(GABA.peak_cf(1:n_GABA,:),...
    numel(GABA.peak_cf(1:n_GABA,:)),1);
all_cf = [all_cf_GABA];
all_cf(all_cf==0) = [];


sumd = zeros(n_runs_km,1); % keep track of sum(distance from centroid)
idx_n = zeros(n_runs_km, length(all_cf));
for n = 1:n_runs_km
    [idx_n(n,:),~,sumd_bands] = kmeans(all_cf,n_bands);
    sumd(n) = mean(sumd_bands);
end
[~,best_cluster] = min(sumd); % Choose kmeans run with lowest sumd
idx = idx_n(best_cluster,:);

median_band = zeros(n_bands,1);
ci_band = zeros(n_bands,2);
figure
hold on
for n = 1:n_bands
    histogram(all_cf(idx==n),'binwidth',1,'facecolor',...
        defaults(mod(n,7)+1,:));
    median_band(n) = median(all_cf(idx==n));
    ci_band(n,1) = prctile(all_cf(idx==n),(100-conf_level)/2);
    ci_band(n,2) = prctile(all_cf(idx==n),(100-conf_level)/2+conf_level);
end

[median_band,sort_ind] = sort(median_band);
ci_band = ci_band(sort_ind,:);
band_freqs = [floor(ci_band(:,1)),ceil(ci_band(:,2))];
xlim([0,max(ci_band(:,2))+5]);
%%
if strcmp(sessions,'GABA')
    aperiodic_signal = GABA.aperiodic_signal;
    fooof_freq = GABA.fooof_freq;
    fooof_spectrum = GABA.fooof_spectrum;
    peak_cf = GABA.peak_cf;
    peak_pow = GABA.peak_pow;

    N = size(GABA.fooof_spectrum,1);
    peak_pow_init= GABA.peak_pow;

    n_sessions = max(session_tags_GABA);
else
    aperiodic_signal = noic.aperiodic_signal;
    fooof_freq = noic.fooof_freq;
    fooof_spectrum = noic.fooof_spectrum;
    peak_cf = noic.peak_cf;
    peak_pow = noic.peak_pow;

    N = size(noic.fooof_spectrum,1);
    peak_pow_init= noic.peak_pow;

    n_sessions = 3;
end



tags = [ones(N/2,1);ones(N/2,1)*2];
tags_bs = [ones(n_bs,1);ones(n_bs,1)*2];
tags_init = tags;

%% Extract peak power and frequency within each band

band_peak_freq = zeros(N,n_bands);
band_peak_pow = zeros(N,n_bands);
peak_pow = 10*peak_pow_init; 
% ^ fooof gives log10(pow), decibels = 10log10(pow)

min_pow = min(min(peak_pow))-1;

for n = 1:N  
    for b = 1:n_bands
        n_freq = peak_cf(n,:);
        n_pow = peak_pow(n,:);
        band_bin = ((n_freq>band_freqs(b,1)).*(n_freq<band_freqs(b,2)))==1;
        if sum(band_bin>=1)
            n_pow(band_bin==0) = min_pow;
            [band_peak_pow(n,b),peak_ind] = max(n_pow);
            band_peak_freq(n,b) = n_freq(peak_ind);
        else
            band_peak_pow(n,b) = 0;
            band_peak_freq(n,b) = NaN;
        end
    end
end

tags = tags_init(1:end-1);

aperiodic = zeros(size(fooof_spectrum));
for n = 1:N
    aperiodic(n,:) = ...
        aperiodic_signal(n,1)-log10(fooof_freq.^(aperiodic_signal(n,2)));
end
aperiodic = 10*aperiodic;


%% Compute statistics

n_freq = size(aperiodic,2);
band_peak_pow_bs = cell(n_bands,2);
band_peak_diff_bs = cell(n_bands,1);
aperiodic_bs = cell(2,2);
aperiodic_diff_bs = cell(2,1);

mean_band_peak_pow_bs = cell(n_bands,2);
mean_band_peak_diff_bs = cell(n_bands,1);
mean_aperiodic_bs = cell(2,2);
mean_aperiodic_diff_bs = cell(2,1);
mean_spectrum_bs = cell(1,2);
mean_spectrum_diff_bs = cell(1,1);

mean_band_diff = zeros(n_bands,1);
ci_band_diff = zeros(n_bands,2);
mean_band_pow = zeros(n_bands,2);
ci_band_pow = zeros(n_bands,2,2);
cf_bs = cell(n_bands,2);
mean_cf_bs = cell(n_bands,2);
mcf = zeros(n_bands,2);
ci_mcf = zeros(n_bands,2,2);
mean_aperiodic_diff = zeros(2,1);
ci_aperiodic_diff = zeros(2,2);
mean_aperiodic = zeros(2,2);
ci_aperiodic = zeros(2,2,2);

for b = 1:n_bands
    for j = 1:2
        ind_start = (j-1)*n_bs*n_samples*n_sessions+1;
        ind_end = min([j*n_bs*n_samples*n_sessions,...
            2*n_bs*n_samples*n_sessions]);
        ind = ind_start:ind_end;
        band_peak_pow_bs{b,j} = band_peak_pow(ind,b);
        cf_bs{b,j} = band_peak_freq(ind,b);
        [~, mean_band_peak_pow_bs{b,j},...
            mean_band_pow(b,j), ci_band_pow(b,j,:)] = ...
            fooof_extract_bs(band_peak_pow(ind,b),n_bs,conf_level);
        
        [~, mcf_bs{b,j},mcf(b,j), ci_mcf(b,j,:)] = ...
            fooof_extract_bs(band_peak_freq(ind,b),n_bs,conf_level);
    end
    band_peak_diff_bs{b} = band_peak_pow_bs{b,2} - ...
        band_peak_pow_bs{b,1};
    [~, mean_band_peak_diff_bs{b},mean_band_diff(b), ...
        ci_band_diff(b,:)] = ...
            fooof_extract_bs(band_peak_diff_bs{b},n_bs,conf_level);
end  



for a = 1:2
    for j = 1:2
        ind_start = (j-1)*n_bs*n_samples*n_sessions+1;
        ind_end = min([j*n_bs*n_samples*n_sessions,...
            2*n_bs*n_samples*n_sessions]);
        ind = ind_start:ind_end;
        aperiodic_bs{a,j} = aperiodic_signal(ind,a);
        [~, mean_aperiodic_bs{a,j},...
            mean_aperiodic(a,j), ci_aperiodic(a,j,:)] = ...
            fooof_extract_bs(aperiodic_signal(ind,a),n_bs,conf_level);
    end
    aperiodic_diff_bs{a} = aperiodic_bs{a,2} - ...
        aperiodic_bs{a,1};
    [~, mean_aperiodic_diff_bs{a},mean_aperiodic_diff(a), ...
        ci_aperiodic_diff(a,:)] = ...
            fooof_extract_bs(aperiodic_diff_bs{a},n_bs,conf_level);
end

%Slightly different structure for spectrum:
mean_spectrum_bs{1} = zeros(n_bs,n_freq);
mean_spectrum_bs{2} = zeros(n_bs,n_freq);
n_freq = size(fooof_spectrum,2);
for a = 1:n_freq
    for j = 1:2
        ind_start = (j-1)*n_bs*n_samples*n_sessions+1;
        ind_end = min([j*n_bs*n_samples*n_sessions,...
            2*n_bs*n_samples*n_sessions]);
        ind = ind_start:ind_end;

        [~, mean_spectrum_bs{j}(:,a)] = ...
            fooof_extract_bs(fooof_spectrum(ind,a),n_bs,conf_level);
    end
    mean_spectrum_diff_bs = mean_spectrum_bs{2} - ...
            mean_spectrum_bs{1};
end


%% Plots

if ~suppress_figs
    % Plot band power + aperiodic components before/after infusion

    ax = [];
    figure
    for b = 1:n_bands
        plot(b-0.2+randn(n_bs,1)*0.05,...
            mean_band_peak_pow_bs{b,1},'.','color',defaults(b,:))
        hold on
        plotConfInterval(b-0.2,mean_band_pow(b,1),...
            ci_band_pow(b,1,1),ci_band_pow(b,1,2),[0,0,0]);
        plot(b+0.2+randn(n_bs,1)*0.05,...
            mean_band_peak_pow_bs{b,2},'.','color',defaults(b,:))
        plotConfInterval(b+0.2,mean_band_pow(b,2),...
            ci_band_pow(b,2,1),ci_band_pow(b,2,2),[0,0,0]);
        ylabel('Mean Power (decibels)')
        xlabels{b} = [num2str(band_freqs(b,1)),'-',num2str(band_freqs(b,2))];

    end
    ax = [ax,gca];
    xlim([0,n_bands+1])
    xticks(1:n_bands);
    xticklabels(xlabels)
    xlabel('Frequency (Hz)')


    aperiodic_labels = {'offset','exponent'};

    for a = 1:2
        figure
        plot(1+randn(n_bs,1)*0.1,...
            mean_aperiodic_bs{a,1},'.','color',defaults(2,:))
        hold on
        plotConfInterval(1,mean_aperiodic(a,1),...
            ci_aperiodic(a,1,1),ci_aperiodic(a,1,2),[0,0,0]);
        plot(2+randn(n_bs,1)*0.1,...
            mean_aperiodic_bs{a,2},'.','color',defaults(2,:))
        plotConfInterval(2,mean_aperiodic(a,2),...
            ci_aperiodic(a,2,1),ci_aperiodic(a,2,2),[0,0,0]);
        ylabel(['Aperiodic ',aperiodic_labels{a}])
        xlim([0,3])
        ax = [ax,gca];


    end
    %  Plot changes in power + aperiodic components

    bands = [1:7];
    figure
    xlabels=cell(1,length(bands));
    ax = [];
    for n = 1:length(bands)
        plot(n+randn(n_bs,1)*0.1,...
            mean_band_peak_diff_bs{bands(n),1},'.','color',defaults(n,:))
        hold on
        plotConfInterval(n,mean_band_diff(bands(n)),...
            ci_band_diff(bands(n),1),ci_band_diff(bands(n),2),[0,0,0]);
        xlabels{n} = [num2str(band_freqs(bands(n),1)),'-',...
            num2str(band_freqs(bands(n),2))];
        ylabel('Mean Change in Power (decibels)')


    end
    ax = [ax,gca];
    xlim([0,length(bands)+1])
    xticks(1:length(bands));
    xticklabels(xlabels)
    %ylim([-5,0.5])
    hold on
    plot([0,length(bands)+1],[0,0],'k--')

    figure
    for n = 1:2
        plot(n+randn(n_bs,1)*0.1,...
            mean_aperiodic_diff_bs{n,1},'.','color',defaults(n,:))
        hold on
        plotConfInterval(n,mean_aperiodic_diff(n),...
            ci_aperiodic_diff(n,1),ci_aperiodic_diff(n,2),[0,0,0]);

        ylabel('Aperiodic parameter')

        ax = [ax,gca];
    end
    xlim([0,3])
    xticks(1:2);
    xticklabels(aperiodic_labels)
    plot([0,3],[0,0],'k--')

    % Plot spectrum

    spectrum = [mean_spectrum_bs{1};mean_spectrum_bs{2}];

    plotSpectrum(10.^(spectrum),fooof_freq',tags_bs,2);
    hold on
    ind_before = tags==1;
    ind_after = tags==2;

    plot(fooof_freq,mean(aperiodic(ind_before,:)),'--','color',...
        defaults(2,:));
    plot(fooof_freq,mean(aperiodic(ind_after,:)),'--','color',...
        defaults(3,:));

    xlim([0,100])
    ylim([-15,20])

    % Remove aperiodic component

    spectrum = [mean_spectrum_bs{1};mean_spectrum_bs{2}];
    tags_bs = [ones(100,1);ones(100,1)*2];
    % Mean aperiodic signal across all samples:
    mean_aperiodic = [mean(aperiodic(ind_before,:));...
        mean(aperiodic(ind_after,:))]
    plotSpectrum_periodic(10.^(spectrum),fooof_freq',tags_bs,2,mean_aperiodic);
    hold on

    xlim([-0,100])
    yl = ylim;
    ylim([-0.5,yl(2)]);
    ylabel('Power above aperiodic component (dB)')
end

%% Output structure
fooof_data.band_freqs = band_freqs;
fooof_data.mean_band_peak_pow_bs = mean_band_peak_pow_bs;
fooof_data.mean_band_pow = mean_band_pow;
fooof_data.ci_band_pow = ci_band_pow;
fooof_data.mean_aperiodic_bs = mean_aperiodic_bs;
fooof_data.mean_aperiodic = mean_aperiodic;
fooof_data.ci_aperiodic = ci_aperiodic;
fooof_data.mean_band_peak_diff_bs = mean_band_peak_diff_bs;
fooof_data.mean_band_diff = mean_band_diff;
fooof_data.ci_band_diff = ci_band_diff;
fooof_data.mean_aperiodic_diff_bs = mean_aperiodic_diff_bs;
fooof_data.mean_aperiodic_diff = mean_aperiodic_diff;
fooof_data.ci_aperiodic_diff = ci_aperiodic_diff;
fooof_data.mean_spectrum_bs = mean_spectrum_bs;
fooof_data.aperiodic = aperiodic;

end