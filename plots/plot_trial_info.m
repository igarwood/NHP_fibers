function [] = plot_trial_info(yl,sample_start,sample_end,...
    reward_start,reward_end,match_start,face_alpha,correct)

    n = correct+2;
    x_fill = [sample_start,sample_end,sample_end,sample_start];
    y_fill = [yl(1),yl(1),yl(2),yl(2)];
    box = patch(x_fill,y_fill,'k','FaceAlpha',face_alpha,'linestyle',...
        'none');
    uistack(box,'bottom') ;
    line = plot(ones(1,100)*match_start,linspace(yl(1),yl(2)),'color',...
        [0.7,0.7,0.7],'linewidth',face_alpha*10);
    uistack(line,'bottom');

    x_fill = [reward_start,reward_end,reward_end,reward_start];
    y_fill = [yl(1),yl(1),yl(2),yl(2)];
    if mod(n,2)==0
        color = [0.4660    0.6740    0.1880];
    else
        color = [0.8500    0.1250    0.0980];
    end
%     color = [0,0,0];
    box = patch(x_fill,y_fill,'k','facecolor',color,...
        'FaceAlpha',face_alpha,'linestyle','none');
    uistack(box,'bottom') ;
    ylim(yl);
end