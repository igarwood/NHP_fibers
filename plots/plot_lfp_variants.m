function raster = plot_lfp_variants(Y, x_time, y_times, ...
    modulation_start, modulation_end, variants, task)
% Plot rasters with the same number of trials for each variant
% For variant(s) with more trials, sample trials from ~roughly the same
% time as the variant with the fewest trials.
modulation_color = [0.65 0.9 0.9210];

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
    
    Y_plot = [Y(:,trials_plot)];
    
    imagesc(x_time,1:size(Y_plot,2),Y_plot')
    hold on
    
    
    for j = 1:length(modulation_start)
        modulation_trials = ...
            get_modulation_trials(y_times(:,trials_plot),...
            modulation_start(j),modulation_end(j));
        x_fill = [x_time(1),x_time(end),x_time(end),x_time(1)];
        y_fill = [modulation_trials(1),modulation_trials(1),...
            modulation_trials(end),modulation_trials(end)];
        box = patch(x_fill,y_fill,modulation_color,'FaceAlpha',...
            0.7,'linestyle','none');
    end
    axis xy
    colormap('gray')
    Y_flat = reshape(Y_plot,1,[]);
    caxis([quantile(Y_flat,0.025),quantile(Y_flat,0.975)])
    
    ax = [ax, gca];
    
end

raster.fig = gcf;
raster.ax = ax;