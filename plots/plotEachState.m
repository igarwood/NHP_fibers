function [ax] = plotEachState(spectrogram, path, spec_freq)

K = max(path);
for k = 1:K
    rows = ceil(K/2);
    cols = 2;
    ax = subplot(rows,cols,k);
    imagesc(1:length(path==k), spec_freq, pow2db(spectrogram(path==k,:,1)'));
    box off
    %set(gca,'clim',[-20 50])
    set(gca,'xticklabel',[])
    axis xy; 
    ylabel('Frequency (Hz)');
    title(strcat('State',num2str(k)));
    ylim([0 100])
    %caxis([-50,0])
    caxis([-20,30])
    cmap = load('white_spec_scale');
    %colormap(gca,'jet')
    colormap(gca,cmap.cmap)
    ax = gca;
end