function [Z,U] = time_rescaling_theorem(lambda,spikes)


spike_ind = find(spikes==1);

Z = zeros(length(spike_ind),1);
if ~isempty(spike_ind)
    Z(1) = sum(lambda(1:spike_ind(1)));


    for n = 2:length(spike_ind)
        Z(n) = sum(lambda(spike_ind(n-1):spike_ind(n)));
    end
end

U = 1-exp(-Z);

end
    