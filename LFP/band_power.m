function [y] = scale_power(spectrogram,spec_freq, freq_bands)

spectrogram_db = pow2db(spectrogram);
H = size(freq_bands,1);
N = size(spectrogram_db,1);
y = zeros(N,H);
y = zeros(N,H);

for h = 1:H
    [~,low] = min(abs(spec_freq-freq_bands(h,1)));
    [~,high] = min(abs(spec_freq-freq_bands(h,2)));
    y(:,h) = mean(spectrogram_db(:,low:high),2)';
end

end