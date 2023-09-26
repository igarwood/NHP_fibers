function fig = plotConfInterval_cont(x,Y,CI_pos,CI_neg,color_ind)
if nargin == 4
    color_ind = 0;
end
K = size(Y,1);


defaults = [ 0    0.4470    0.7410 
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

for k = 1:K
    y = Y(k,:);
    x_fill =[x, fliplr(x)];
    y_fill =[CI_neg(k,:), fliplr(CI_pos(k,:))];

    if k < 8
        plot(x,y,'Color',defaults(k+color_ind,:));
        %semilogy(sfreqs,(spectrum),'Color',defaults(i,:));
        hold on
        fill(x_fill,y_fill,defaults(k+color_ind,:),'FaceAlpha',0.5,'EdgeColor',...
            defaults(k+color_ind,:),'linestyle','none');
    elseif k < 15
        plot(x,y,'Color',defaults(k-7,:));
        hold on
        fill(x_fill,y_fill,defaults(k-7,:),'FaceAlpha',0.5,'EdgeColor',...
            defaults(k-7,:),'linestyle','none');
    else
        plot(x,y,'Color',defaults(k-14,:));
        hold on
        fill(x_fill,y_fill,defaults(k-14,:),'FaceAlpha',0.5,'EdgeColor',...
            defaults(k-14,:),'linestyle','none');
    end
    
end