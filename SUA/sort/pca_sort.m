function spikes_cluster = pca_sort(spike_dur, unit_channel, num_spikes,...
   nchan)

% Compute PCs from smoothed spike waveforms; this approach
% helps to capture slower dynamics in the spike waveforms
ds_fact = 3;
ds = round(spike_dur/ds_fact);
unit_channel_ds = cell(4,1);
for c = 1:nchan
    for d = 1:ds
        unit_channel_ds{c}(:,d) = mean(unit_channel{c}(:,ds_fact*(d-1)+1:ds_fact*d),2);
    end
end

X = [];
for n = 1:nchan
    X = [X,unit_channel_ds{n}];
end



[eig_vec, eig_val] = eig(cov(X));
[eig_val,sort_ind] = sort(diag(eig_val),'descend');
eig_vec = eig_vec(:,sort_ind);
scores = X*eig_vec;

% Additionally compute PCs from the first 20 samples of raw spike waveforms
ds_fact = 1;
ds = round(20/ds_fact);
unit_channel_ds = cell(4,1);
for c = 1:nchan
    for d = (1:ds)+6
        unit_channel_ds{c}(:,d-6) = mean(unit_channel{c}(:,ds_fact*(d-1)+1:ds_fact*d),2);
    end
end

X = [];
for n = 1:nchan
    X = [X,unit_channel_ds{n}];
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