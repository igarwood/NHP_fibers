function band_freqs = fooof_bands(n_bands, conf_level,n_runs_km)
% Function to extract frequency bands from PMC GABA and saline data

% Parameters from FOOOF analysis:
max_peaks = 6; 

% Post-processing parameters:
if nargin < 1
    n_bands = 7;
end
if nargin < 2;
    conf_level = 95;
end
if nargin < 3
    n_runs_km = 10;
end

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
folder = [info.save_folder,'/data_for_fooof/'];
fileroot = 'fooof_output';

%% Load data

GABA = load([folder,fileroot,'_GABA_',num2str(max_peaks),...
    '.mat']);
saline = load([folder,fileroot,'_saline_',num2str(max_peaks),...
    '.mat']);

spec_sample_GABA = readmatrix([folder,'spec_sample_GABA.csv']);
spec_sample_saline = readmatrix([folder,'spec_sample_saline.csv']);

spec_freq = readmatrix([folder,'freqs.csv']);
load([folder,'session_tags'])

n_GABA = size(GABA.fooof_spectrum,1);
n_saline = size(saline.fooof_spectrum,1);

%% Determine frequency bands of interest:

all_cf_saline = reshape(saline.peak_cf(1:n_saline,:),...
    numel(saline.peak_cf(1:n_saline,:)),1);
all_cf_GABA = reshape(GABA.peak_cf(1:n_GABA,:),...
    numel(GABA.peak_cf(1:n_GABA,:)),1);
all_cf = [all_cf_saline;all_cf_GABA];
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

end