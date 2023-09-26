function raster = plot_raster_variants(Y, x_time, y_times, ...
    modulation_start, modulation_end, variants, task)
% Plot rasters with the same number of trials for each variant
% For variant(s) with more trials, sample trials from ~roughly the same
% time as the variant with the fewest trials.


n_variants = length(variants);
for n = 1:n_variants
    n_trials_variant(n)=sum(task==variants(n));
end
[min_n_trials, min_variant] = min(n_trials_variant);
min_variant_ind = find(task==variants(min_variant));

ax = [];
figure
for n = 1:n_variants
    subplot(1,n_variants,n)
    
    trials_plot = zeros(1,min_n_trials);
    variant_ind = find(task == variants(n));
    for m = 1:min_n_trials
        [~,min_ind] = min(abs(variant_ind-min_variant_ind(m)));
        trials_plot(m) = variant_ind(min_ind);
        variant_ind(min_ind) = [];
    end
    
    trials_plot = sort(trials_plot);
    
    modulation_trials = ...
        get_modulation_trials(y_times(:,trials_plot),modulation_start,...
        modulation_end);
    if isempty(modulation_trials)
        raster_var = plot_raster(Y(:,trials_plot), x_time);
    else
        raster_var = plot_raster(Y(:,trials_plot), x_time, [], ...
            modulation_trials(1), modulation_trials(2));
    end
    ax = [ax, raster_var.ax];
    
end

raster.fig = gcf;
raster.ax = ax;