function raster = plot_all_rasters(Y, x_time, y_times, modulation_start,...
    modulation_end, task_variant, variants, variant_names, task,time_yn)
% time_yn = 1, y-axis of relevant rasters will be time instead of trials

if time_yn == 1
    y_time = mean(y_times)/60; % in minutes
else
    y_time = [];
end

load('trial_info.mat')

raster = cell(1,2);
%% Plot raster for all trials

modulation_trials = get_modulation_trials(y_times,...
    modulation_start,modulation_end);

figure
raster{1} = plot_raster(Y, x_time, y_time, modulation_trials(1), ...
    modulation_trials(2));
title('All trials')
yl = ylim;
plot_trial_info(yl,sample_start,sample_end,...
    reward_start,reward_end,match_start,0.2,1)
xlabel('Time since cue onset (seconds)')
ylabel('Trial #')

%% Plot raster for all trials from each variant
n_variants = length(variants);

for n = 1:n_variants
    figure
    trials_plot = task==variants(n);
    modulation_trials = ...
        get_modulation_trials(y_times(:,trials_plot),modulation_start,...
        modulation_end);
    raster{1+n} = plot_raster(Y(:,trials_plot), x_time,...
        y_time(trials_plot), modulation_trials(1), modulation_trials(2));
    title(['All ', variant_names{n},' trials'])
    yl = ylim;
    
    if strcmp(task_variant,'rewarded')
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,variants(n));
    else
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1);
    end
    xlabel('Time since cue onset (seconds)')
    ylabel('Trial #')
end
%% Plot raster for matched trials from each variant 
% Plot the same number of trials for each variant
raster{2+n_variants} = plot_raster_variants(Y, x_time, y_times, ...
    modulation_start, modulation_end, variants, task);

for n = 1:n_variants
    subplot(1,n_variants,n)
    hold on
    yl = ylim;
    title([variant_names{n},' trials'])
    if strcmp(task_variant,'rewarded')
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,variants(n));
    else
        plot_trial_info(yl,sample_start,sample_end,...
            reward_start,reward_end,match_start,0.2,1);
    end
    xlabel('Time since cue onset (seconds)')
    if n == 1
        ylabel('Trial #')
    end
end

end