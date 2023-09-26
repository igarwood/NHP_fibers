function figs = plot_spike_waves(num_spikes,spikes,samp_spikes,channels,...
    spike_time)

nchan = length(channels);

for spike_ind = 1:num_spikes
    figure('Name',['Unit ',num2str(spike_ind)]);
    set(gcf,'renderer','Painters')
    ax = [];
    for c = 1:nchan
        subplot(length(channels),1,c)
        spike = spikes{c,spike_ind};

        if length(spike)>500
            samp_spikes= randsample(length(spike),500);
        else
            samp_spikes= 1:size(spike,1);
        end

        mean_spike = mean(spike);
        plot(spike_time*1000,spikes{c,spike_ind}(samp_spikes,:),'color',...
            [0,0,0,0.02])%[0, 0.4470, 0.7410 0.01])
        hold on
        plot(spike_time*1000,mean_spike,'color',...
            [0,0,0],'linewidth',2)
        ax = [ax,gca];
        xlim([0,2])
        xlabel('Time (ms)')
        ylabel('Amplitude (microvolts)')
    end
    linkaxes(ax,'xy')
    figs(spike_ind) = gcf;
end



end