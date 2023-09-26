function [spikes,spikes_locs,samp_spikes] = assign_spikes(...
    spikes_cluster,unit_channel,unit_locs_channel,N,nchan)
    
num_spikes = max(spikes_cluster);

spikes = cell(nchan,num_spikes);
spikes_locs = cell(1,num_spikes);
samp_spikes = cell(nchan,num_spikes);
for j = 1:num_spikes
    for k = 1:nchan
        spikes{k,j} = unit_channel{k}(spikes_cluster==j,:);
        if length(spikes{k,j})>500
            samp_spikes{k,j} = randsample(length(spikes{k,j}),500);
        else
            samp_spikes{k,j} = 1:size(spikes{k,j},1);
        end
    end
    spikes_locs{j} = unit_locs_channel{1}(spikes_cluster==j);
end

end