function fig = plot_spike_ts(NS6,spikes_locs,fs,spike_dur,win,channels,...
    chan_plot,unit_plot,N_n,w)

figure
set(gcf,'renderer','Painters')

time = 0:1/fs:(N_n-1)/fs;
ind = floor(spike_dur*0.2);

time_plot = (time + ((win-1)*N_n)/fs);



chan_ind = zeros(size(chan_plot));
for c = 1:length(chan_plot)
    chan_ind(c) = find(channels == chan_plot(c));
end
ax2 = [];
data_sum =zeros(1,N_n);


colors = [[0, 0.4470, 0.7410];...
    [0.8500, 0.3250, 0.0980];...
    [0.9290, 0.6940, 0.1250];...
    [0.4940, 0.1840, 0.5560];...
    [0.4660, 0.6740, 0.1880];...
    [0.3010, 0.7450, 0.9330];...
    [0.6350, 0.0780, 0.1840];...
    [0, 0.4470, 0.7410]];	

for c = 1:length(chan_plot)
    subplot(length(chan_plot),1,c)
    chan = chan_plot(c);
    data = NS6.Data(chan,((win-1)*N_n+1):(win*N_n))/4;
    
    spike_ts_chan(1,:) = filtfilt(w,1,double(data));
    %spike_ts_chan(1,:) = double(data);
    data_sum = spike_ts_chan(1,:)+data_sum;

    plot(time_plot,spike_ts_chan(1,:),'k');
    hold on
    ax2 = [ax2, gca];
    
    
    subplot(length(chan_plot),1,c)
    for j = 1:length(unit_plot)
        u = unit_plot(j);
        ind_2 = (spikes_locs{u}>=((win-1)*N_n+1)).*((spikes_locs{u}<(win*N_n)));
        ind_2 = ind_2==1;
        plot(time_plot(spikes_locs{u}(ind_2)-(win-1)*N_n),...
            spike_ts_chan(1,spikes_locs{u}(ind_2)-(win-1)*N_n),'x','color',colors(u,:));
        title(['Electrode ',num2str(chan)])
    end
    
    xlabel('Time (seconds)')
    ylabel('Amplitude (microvolts)')
    %ylim([-120,50])
end
linkaxes(ax2,'x');
fig = gcf;