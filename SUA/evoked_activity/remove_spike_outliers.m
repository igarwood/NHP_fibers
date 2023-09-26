function trial_outliers = remove_spike_outliers(Y,n_per_seg)

% Remove trials where the firing rate is significantly different from a
% group of neighboring trials

n_trials = size(Y,2);
n_trial_groups = 10;
n_trials_per_group = floor(n_trials/n_trial_groups);
Y_trial_outliers = zeros(1,n_trials);
for n = 1:(n_trial_groups+1)
    ind = ((n-1)*n_trials_per_group+1):min(n_trials,n*n_trials_per_group);

    Y_reshape = reshape(Y(:,ind,:)',length(ind),[],n_per_seg);
    Y_means = mean(Y_reshape,3);
    Y_means_flattened = reshape(Y_means,1,[]);

    outliers = quantile(Y_means_flattened,0.999);

    if n_per_seg == 2
        Y_means_outliers = (Y_means>0.5);
    else
        Y_means_outliers = (Y_means>outliers);
    end
   
    Y_trial_outliers(ind)= sum(Y_means_outliers,2)>0;

end
trial_outliers = Y_trial_outliers==1;

    
end