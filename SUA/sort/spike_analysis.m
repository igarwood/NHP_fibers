%% Mix of automated and manual spike sorting; run section-by-section
% Also includes some basic analysis + plotting
session = 'pmc_gaba3';

%% 1. Parameters
num_spikes_init = 3;
sort_yn = 0;
save_sort_yn = 1;
noise_thresh = 300;

%% 2. Extract session information:
[session_info] = setup(session);
folder = session_info.folder;
save_folder = session_info.save_folder;
fileroot = session_info.fileroot;
pharmacology = session_info.pharmacology;
fs = session_info.fs_spike;

spike_folder = [save_folder,'spike_data/'];
spike_file = [spike_folder,'spike_data_stereo_',session];
infofile = [folder,'chanmap_',session,'.m'];

[elecs,lfp_elec] = get_elecs(session);
n_elecs = length(elecs);

infusion_start = pharmacology.infusion_start;
infusion_end = pharmacology.infusion_end;
infusion_rate = pharmacology.infusion_rate;

%% 3. Load data and set time variables
load(spike_file);

N_n = 10000000;
N = max(unit_locs{1,end})+N_n*(length(unit_locs)-1);
time = 0:1/fs:(N-1)/fs;

nchan = size(units,1);
channels = [1:nchan];

N_win = length(unit_locs);
spike_ts_chan = zeros(nchan,N_n);
spike_dur = round(0.0025*fs);  
min_spike_dist = round(0.0005*fs);

s_time = time;

spike_time = time(1:spike_dur);


%% 4. Combine windows from spike extraction
unit_channel = cell(nchan,1);
unit_locs_channel = cell(1,1);

for n = 1:N_win
    for c = 1:nchan
        unit_channel{c} = [unit_channel{c}; units{c,n}];
    end
    unit_locs_channel{1} = [unit_locs_channel{1},...
        unit_locs{1,n}+(n-1)*N_n];
end


%% 5. Remove artifact:
for c = 1:nchan
    spike_keep = min(unit_channel{c},[],2) > -1*noise_thresh;

    for i = 1:length(unit_channel)
        unit_channel{i} = unit_channel{i}(spike_keep,:);
    end
    unit_locs_channel{1} = unit_locs_channel{1}(spike_keep);
end

for c = 1:nchan
    spike_keep = max(unit_channel{c},[],2) < noise_thresh;

    for i = 1:length(unit_channel)
        unit_channel{i} = unit_channel{i}(spike_keep,:);
    end
    unit_locs_channel{1} = unit_locs_channel{1}(spike_keep);
end

%% 6. Sorting based on PCA
num_spikes = num_spikes_init;
if sort_yn
    spikes_cluster = pca_sort(spike_dur, unit_channel, num_spikes,nchan);
else
    load([spike_folder,'spikes_cluster_',session]);
end

%% 6.1 [Manual Step] Recombine units 
units_combine = [3,5];
units_combine = sort(units_combine);
for n = 2:length(units_combine)
    spikes_cluster_orig = spikes_cluster;
    orig_ind = find(spikes_cluster_orig==units_combine(n));
    spikes_cluster(orig_ind) = units_combine(1);
    for m = (units_combine(n)+1):num_spikes
        orig_ind = find(spikes_cluster_orig==m);
        spikes_cluster(orig_ind) = m-1;
    end
    num_spikes = num_spikes-1;
end

%% 6.2 [Manual Step] Resort with a subset of channels 
% Can be useful if unit is predominantly in a subset of channels)

unit_resort = 1;
chan_subset=[4];
orig_ind = find(spikes_cluster==unit_resort);

spikes_cluster_resort = pca_sort_subset(spike_dur, unit_channel,...
    2,chan_subset,unit_resort,spikes_cluster);
num_spikes = num_spikes+1;

spikes_cluster(orig_ind(spikes_cluster_resort==1)) = unit_resort;
spikes_cluster(orig_ind(spikes_cluster_resort==2)) = num_spikes;

%% 7. Assign + plot spike waveforms
N = N_n*N_win;
num_spikes = max(spikes_cluster);
[spikes,spikes_locs,samp_spikes] = ...
    assign_spikes(spikes_cluster,unit_channel,unit_locs_channel,N,nchan);
figs = plot_spike_waves(num_spikes,spikes,samp_spikes,channels,...
        spike_time);
    
if save_sort_yn
    save([spike_folder,'spikes_cluster_',session],'spikes_cluster','spikes_locs');
end

%% 8. Plot mean spike rates
figs = [figs,plot_spike_rate(N, num_spikes,N_win,spikes,spikes_locs,...
    time, infusion_start,infusion_end)];

%% 9. Plot spike waveforms before and after infusion
figs = [figs,plot_spikes_before_after(num_spikes,spikes,spikes_locs,...
    channels,spike_time,infusion_start,infusion_end)];

%% 10. Plot spike waveforms during infusion
figs = [figs,plot_spikes_during(num_spikes,spikes,spikes_locs,...
    channels,spike_time,start_inject,end_inject)];

