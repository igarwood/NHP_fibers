function spikes_cluster = pca_sort_subset(spike_dur, unit_channel,...
    num_spikes,chan_subset,unit_resort,spikes_cluster_init)

nchan = length(chan_subset);
if nargin == 6
    ind = spikes_cluster_init==unit_resort;
else
    ind = 1:length(unit_channel{1});
end
ds_fact = 3;
ds = round(spike_dur/ds_fact);
unit_channel_ds = cell(4,1);

X = [];
for c = 1:nchan
    chan = chan_subset(c);
    for d = 1:ds
        unit_channel_ds{c}(:,d) = ...
            mean(unit_channel{chan}(ind,ds_fact*(d-1)+1:ds_fact*d),2);
    end
    X = [X,unit_channel_ds{c}];
end


[eig_vec, eig_val] = eig(cov(X));
[eig_val,sort_ind] = sort(diag(eig_val),'descend');
eig_vec = eig_vec(:,sort_ind);
scores = X*eig_vec;

ds_fact = 1;
ds = round(20/ds_fact);
unit_channel_ds = cell(4,1);
for c = 1:nchan
    chan = chan_subset(c);
    for d = (1:ds)+6
        unit_channel_ds{c}(:,d-6) = ...
            mean(unit_channel{chan}(ind,ds_fact*(d-1)+1:ds_fact*d),2);
    end
    X = [X,unit_channel_ds{c}];
end


[eig_vec, eig_val] = eig(cov(X));
[eig_val,sort_ind] = sort(diag(eig_val),'descend');
eig_vec = eig_vec(:,sort_ind);
scores_short = X*eig_vec;

obs = zeros(size(scores,1),10);
obs(:,1:2:10) = scores(:,1:5);
obs(:,2:2:10) = scores_short(:,1:5);


spikes_cluster = kmeans(obs,num_spikes,'maxiter',1000);

    

end
