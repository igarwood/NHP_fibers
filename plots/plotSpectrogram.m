function [ax] = plotSpectrogram(spectrogram, time, spec_freq,tunit)


if nargin<4 tunit = 'minutes'; end

imagesc(time, spec_freq, pow2db(spectrogram(:,:,1)'));
ax = gca;
set(gca, 'Ticklength', [0 0])

imagesc(time, spec_freq, pow2db(spectrogram(:,:,1)'));
ax = gca;
set(gca, 'Ticklength', [0 0])

box off
%set(gca,'clim',[-20 50])
axis xy; 
ylabel('Frequency (Hz)');
% c = colorbar;
% ylabel(c,'Power (dB)');
ylim([0 305])
caxis([-20,30])
set(gcf,'renderer','Painters')

cmap = load('white_spec_scale');
%colormap(gca,'jet')
colormap(gca,cmap.cmap)
xlim([min(time),max(time)]);

if strcmp(tunit, 'minutes')
    xlabel('Time (min)')
else
    xlabel('Time (s)')
end

end