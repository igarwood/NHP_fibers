function figs = plot_spikes_during(num_spikes,spikes,spikes_locs,...
    channels,spike_time,start_inject,end_inject)

channels = channels-min(channels)+1;
for spike_ind = 1:num_spikes
    spike = zeros(size(spikes{1,spike_ind},1),length(spike_time),...
        length(channels));
    for c = 1:length(channels)
        spike(:,:,c) = spikes{c,spike_ind};
    end
    
    spike_locs = spikes_locs{spike_ind};
    
    N = max(spike_locs);%NS6.MetaTags.DataPoints;
    fs = 30*10^3;
    time = 0:1/fs:(N-1)/fs;
    spike_times = time(spike_locs);
    
    spike_ind_during = ((spike_times>start_inject(1)).*(spike_times<end_inject(1)))==1;

    
    spike_locs_during = spike_locs(spike_ind_during);
    
    spike_during = spike(spike_ind_during,:,:);
    
    
    try
        samp_spike_during = randsample(length(spike_during),200);
    catch 
        samp_spike_during = 1:size(spike_during,1);
    end

 
    
    mean_spike_during = mean(spike_during);
    
    figure('Name',['Unit ',num2str(spike_ind)]);
    set(gcf,'renderer','Painters')
    figs(spike_ind) = gcf;
    ax = [];
    set(gcf,'renderer','Painters')
    for c = 1:length(channels)
        chan = channels(c);
        subplot(length(channels),1,chan)
        try
            plot(spike_time*1000,spike_during(samp_spike_during,:,c),'color',...
                [0,0,0,0.02])%[0, 0.4470, 0.7410 0.01])
            hold on
            plot(spike_time*1000,mean_spike_during(:,:,c),'color',...
                [0,0,0],'linewidth',2)
        end
        xlim([0,2.5])
        if chan == 1
            title('Before injection')
        elseif chan == length(channels)
            xlabel('Time (ms)')
        end
        ylabel('Amplitude (microvolts)')
        ax = [ax,gca];


        linkaxes(ax,'y')
    end
end