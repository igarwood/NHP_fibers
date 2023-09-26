function [all_samples, mean_samples, median_bs, ci_bs] = ...
    fooof_extract_bs(samples,n_bs,conf_level)

% 
N_samples = size(samples,1);
all_samples = zeros(floor(N_samples/n_bs),n_bs);
for n = 1:n_bs
    bs_ind= randsample(N_samples,floor(N_samples/n_bs),'true');
    all_samples(:,n) = samples(bs_ind,:);
end


%all_samples = reshape(samples,[],n_bs);
mean_samples = mean(all_samples,'omitnan');

median_bs = median(mean_samples);
ci_bs(1) = prctile(mean_samples,(100-conf_level)/2);
ci_bs(2) = prctile(mean_samples,(100-conf_level)/2+conf_level);

end