fooof_data_GABA = fooof_analysis_put('GABA',1);
fooof_data_extra = fooof_analysis_put('noic',1);
%%

% Default colors for plots:
defaults = [0.6350    0.0780    0.1840;...
        0.0000    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330];

ax = [];
n_bs = 100;

%  Plot changes in power + aperiodic components

figure
bands = [1:7];
xlabels={'GABA','no inf'};
ax = [];
for n = 1:length(bands)
    subplot(1,9,n)
    plot(1+randn(n_bs,1)*0.1,...
        fooof_data_GABA.mean_band_peak_diff_bs{bands(n),1},'.','color',defaults(n,:))
    hold on
    plotConfInterval(1,fooof_data_GABA.mean_band_diff(bands(n)),...
        fooof_data_GABA.ci_band_diff(bands(n),1),fooof_data_GABA.ci_band_diff(bands(n),2),[0,0,0]);
    
    plot(2+randn(n_bs,1)*0.1,...
        fooof_data_extra.mean_band_peak_diff_bs{bands(n),1},'.','color',defaults(n,:))
    plotConfInterval(2,fooof_data_extra.mean_band_diff(bands(n)),...
        fooof_data_extra.ci_band_diff(bands(n),1),fooof_data_extra.ci_band_diff(bands(n),2),[0,0,0]);
    
    ylabel('Change in Mean Power (decibels)')
    xlim([0,3])
    xticks(1:3);
    xticklabels(xlabels)
    ax = [ax,gca];
    
    plot([0,3],[0,0],'k--')

end

linkaxes(ax,'y');
ylabels = {'offset','exponent'}
ax = [];
for n = 1:2
    subplot(1,9,n+7)
    plot(1+randn(n_bs,1)*0.1,...
        fooof_data_GABA.mean_aperiodic_diff_bs{n,1},'.','color',defaults(n,:))
    hold on
    plotConfInterval(1,fooof_data_GABA.mean_aperiodic_diff(n),...
        fooof_data_GABA.ci_aperiodic_diff(n,1),fooof_data_GABA.ci_aperiodic_diff(n,2),[0,0,0]);
    
    plot(2+randn(n_bs,1)*0.1,...
        fooof_data_extra.mean_aperiodic_diff_bs{n,1},'.','color',defaults(n,:))
    plotConfInterval(2,fooof_data_extra.mean_aperiodic_diff(n),...
        fooof_data_extra.ci_aperiodic_diff(n,1),fooof_data_extra.ci_aperiodic_diff(n,2),[0,0,0]);
    

    ylabel(ylabels{n})
    xlim([0,3])
    xticks(1:3);
    xticklabels(xlabels)
    ax = [ax,gca];
    
    plot([0,3],[0,0],'k--')
end
linkaxes(ax,'y');