function fig = plotSpectrum(spectrogram,spec_freq,path,K)
figure

% Default colors:
defaults = [
        0.0000    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.6350    0.0780    0.1840];
%     
% freq_keep = ((spec_freq <59.6) + (spec_freq > 60.4)) == 1;
% spec_freq = spec_freq(freq_keep);
% spectrogram = spectrogram(:,freq_keep);
for k = 1:K
    
    spec = pow2db(spectrogram(path==k,:));
    spectrum = mean(spec);
    CI = zeros(2,length(spectrum));
    for n = 1:length(spectrum)
        CI(1,n) = prctile(spec(:,n),2.5);
        CI(2,n) = prctile(spec(:,n),97.5);
    end
%     CI(1,:) = spectrum - 1.96*std(spec);%/sqrt(size(spec,1));
%     CI(2,:) = spectrum + 1.96*std(spec);%/sqrt(size(spec,1));
    x =[spec_freq, fliplr(spec_freq)];
    y =[CI(1,:), fliplr(CI(2,:))];
    
    if k < 8
        plot(spec_freq,spectrum,'Color',defaults(k,:));
        %semilogy(sfreqs,(spectrum),'Color',defaults(i,:));
        hold on
        fill(x,y,defaults(k,:),'FaceAlpha',0.3,'EdgeColor',...
            defaults(k,:));
    elseif k < 15
        plot(spec_freq,spectrum,'Color',defaults(k-7,:));
        hold on
        fill(x,y,defaults(k-7,:),'FaceAlpha',0.3,'EdgeColor',...
            defaults(k-7,:));
    else
        plot(spec_freq,spectrum,'Color',defaults(k-14,:));hold on
        fill(x,y,defaults(k-14,:),'FaceAlpha',0.3,'EdgeColor',...
            defaults(k-14,:));
    end
    
    xlim([0.5,500]);

end

box off
set(gca, 'Ticklength', [0 0])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
ax = gca;
ax.FontSize = 15;
fig = gcf;


end