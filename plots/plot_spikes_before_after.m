function figs = plot_spikes_before_after(num_spikes,spikes,spikes_locs,...
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
    
    spike_ind_pre = (spike_times<start_inject(1)).*...
        (spike_times>(start_inject(1)-300))==1; 
    spike_ind_post = spike_times>(spike_times(end)-300);%spike_times>end_inject(end);
    
    spike_locs_pre = spike_locs(spike_ind_pre);
    spike_locs_post = spike_locs(spike_ind_post);
    
    spike_pre = spike(spike_ind_pre,:,:);
    spike_post = spike(spike_ind_post,:,:);
    
    
    try
        samp_spike_pre = randsample(length(spike_pre),200);
    catch 
        samp_spike_pre = 1:size(spike_pre,1);
    end

    try
        samp_spike_post = randsample(length(spike_post),200);
    catch 
        samp_spike_post = 1:size(spike_post,1);
    end
%     
    
    mean_spike_pre = mean(spike_pre);
    mean_spike_post = mean(spike_post);
    
    figure('Name',['Unit ',num2str(spike_ind)]);
    set(gcf,'renderer','Painters')
    figs(spike_ind) = gcf;
    ax = [];
    set(gcf,'renderer','Painters')
    for c = 1:length(channels)
        chan = channels(c);
        subplot(length(channels),2,(chan-1)*2+1)
        try
            plot(spike_time*1000,spike_pre(samp_spike_pre,:,c),'color',...
                [0,0,0,0.02])%[0, 0.4470, 0.7410 0.01])
            hold on
            plot(spike_time*1000,mean_spike_pre(:,:,c),'color',...
                [0,0,0],'linewidth',2)
        end
        xlim([0,1.5])
        %ylim([-60,40])
        if chan == 1
            title('Before injection')
        elseif chan == length(channels)
            xlabel('Time (ms)')
        end
        ylabel('Amplitude (microvolts)')
        ax = [ax,gca];


        subplot(length(channels),2,(chan-1)*2+2)
        try
            plot(spike_time*1000,spike_post(samp_spike_post,:,c),'color',...
                [0,0,0,0.02])%[0, 0.4470, 0.7410 0.01])
            hold on
            plot(spike_time*1000,mean_spike_post(:,:,c),'color',...
                [0,0,0],'linewidth',2)
        end
        xlim([0,1.5])
        ylim([-60,40])
        if chan == 1
            title('After injection')
        elseif chan == length(channels)
            xlabel('Time (ms)')
        end
        ax = [ax,gca];
        linkaxes(ax,'y')
    end
end