function raster = plot_raster(Y, x_time, y_time, modulation_start, ...
    modulation_end)
% y-axis can either be trials or trial time; if y_time is empty, use trials
% modulation_start/modulation_end are the start/end trials of modulation
modulation_color = [0.2,0.8,1];
face_alpha = 0.2;


n_trials = size(Y,2);
if nargin == 2
    y_time = [];
    modulation_start = [];
    modulation_end = [];
elseif nargin == 3
    modulation_start = [];
    modulation_end = [];
end

for k = 1:n_trials
    if isempty(y_time)
        plot(x_time(Y(:,k)==1), ones(sum(Y(:,k)==1),1)*k,'.k')
    else
        plot(x_time(Y(:,k)==1),ones(sum(Y(:,k)==1),1)*y_time(k),'.k')
    end
    hold on
end
if isempty(y_time)
    yl = [0,k];
else
    yl = [min(y_time),max(y_time)];
end
xl = [min(x_time),max(x_time)];

ylim(yl)
xlim(xl)

for n = 1:length(modulation_start)
    x_fill = [xl(1), xl(2), xl(2), xl(1)];
    
    y_fill = [modulation_start(n), modulation_start(n),...
        modulation_end(n), modulation_end(n)];
    if ~isempty(y_time)
        y_fill = y_time(y_fill);
    end
    box = patch(x_fill,y_fill,'k','facecolor',modulation_color,...
        'facealpha',face_alpha,'linestyle','none');
    uistack(box,'bottom') ;
end

raster.fig = gcf;
raster.ax = gca;

end

