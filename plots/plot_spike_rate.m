function figs = plot_spike_rate(N, num_spikes,N_win,spikes,spikes_locs,...
    time,start_inject,end_inject)

figure
t=tiledlayout(num_spikes,1, 'Padding', 'none', 'TileSpacing', 'compact');
ax = [];
fs = 1/(time(2)-time(1));
for spike_ind = 1:num_spikes
    time_scale = 60;
    spike_locs = spikes_locs{spike_ind};
    spike_bin= zeros(size(time));
    spike_bin(spike_locs) = 1;

    win_time = 60; %seconds
    win = fs*win_time;

    spike_rate = movmean(spike_bin,win)*fs;
    
    nexttile(t)
    if spike_ind > 1
        set(gca,'XTickLabel',[]);
        set(gca,'XLabel',[]);
    else
        xlabel('Time (min)')
    end
    
    
    %subplot(num_spikes,1,spike_ind)
    plot(time(1:fs:end)/time_scale,spike_rate(1:fs:end))
    hold on
    inject_times = [(start_inject);end_inject]/time_scale;
    yl = ylim;
    plot_loc = yl(2);
    plot(inject_times,ones(size(inject_times))*plot_loc,'k',...
        'linewidth',5)
    xlim([0,max(time)/time_scale])
    set(gca,'YLabel',[]);
    %ylabel('Spike rate (spikes/second)')
    ax = [ax,gca];
    
    title(['Unit ',num2str(spike_ind)])
    
end
set(gcf,'renderer','Painters')
figs(spike_ind) = gcf;
linkaxes(ax,'x')


end