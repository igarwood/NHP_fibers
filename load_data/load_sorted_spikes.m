function [spike_data] = load_sorted_spikes(session,suppress_figs);

%% 1. Parameters
noise_thresh = 300;

%% 2. Extract session information:
[session_info] = setup(session);
folder = session_info.folder;
fileroot = session_info.fileroot;
fs = session_info.fs_spike;
spike_path = [fileparts(matlab.desktop.editor.getActiveFilename)];
slashes = strfind(spike_path,'/');
spike_path = spike_path(1:slashes(end-1));

file = [folder,fileroot,'.ns6'];
spike_folder = [spike_path,'spike_data/'];
spike_file = [spike_folder,'spike_data_stereo_',session];
infofile = [folder,'chanmap_',session,'.m'];

[elecs,lfp_elec] = get_elecs(session);
n_elecs = length(elecs);


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

%% Load spike clusters
load([spike_folder,'spikes_cluster_',session],'spikes_cluster');

%% 7. Assign + plot spike waveforms
N = N_n*N_win;
num_spikes = max(spikes_cluster);
[spikes,spike_locs,samp_spikes] = ...
    assign_spikes(spikes_cluster,unit_channel,unit_locs_channel,N,nchan);
if ~suppress_figs
    figs = plot_spike_waves(num_spikes,spikes,samp_spikes,channels,...
            spike_time);
end
    
%% Output data
spike_data.spikes = spikes;
spike_data.spike_locs = spike_locs;
spike_data.s_time = s_time;

end