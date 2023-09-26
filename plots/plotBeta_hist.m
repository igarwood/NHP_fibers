function fig = plotBeta_hist(Y,path,beta_a, beta_b)
% Plot frequency and state specific frequency bands
% Beta distribution: p(y) = (y^(a-1)*(1-y)^(b-1))/Beta(a,b)
% Input: for H frequency bands and K states,
%        Y = observations (H bands)
%        path = state path where each element is {1,2,...,k}
%        beta_a: KxH matrix of a 
%        beta_b: KxH matrix of b

% Default colors:
defaults = [0.6953    0.1328    0.1328;...
    0.0000    0.0000    0.7500;...
    0.8516    0.6445    0.1250;...
    0.0550    0.6053    0.5441;...
    0.4940    0.1840    0.5560];
defaultLine = max(zeros(size(defaults)),defaults-0.2);

for k = 1:K
    fig(k) = figure('Name',strcat('State',num2str(k)));
end
for j = 1:H
    for k = 1:K
        set(0, 'CurrentFigure', fig(k))
        set(gcf,'renderer','Painters')
        subplot(1,5,j)
        y = Y(path==k,j);
        y= y*(1-4E-12)+2E-12;
        pd = makedist('Beta','a',beta_a(k,j),'b',beta_b(k,j));
        x = linspace(0,1,10000);

        [f,xhist] = ecdf(y);
        ecdfhist(f,xhist,0.015:0.03:0.985);
        h = findobj(gca,'Type','patch');
        h.FaceColor = defaults(j,:);
        h.FaceAlpha = 0.5;
        
        set(gca,'TickLength',[0 0])
        hold on
        ax = gca;
        ax.YColor = 'none';
        c = ax.Color;
        xlim([0,1])
        % Adjust y lims as needed here:
        ylim([0,20])
        box off
        if j == 3
            xlabel('Normalized Power')
        end
        
        Y = pdf(pd,x);
        plot(x,Y,'k-','LineWidth',3,'Color',defaultLine(j,:));
        set(gca,'TickLength',[0 0])
        set(gca,'fontsize',16)
        hold on
        ax = gca;
        ylim([0,8])
        ax.YColor = 'k';
        ax = gca;
        if j == 1
            ylabel('f(x)')
        else
            set(gca,'yticklabel',[])
        end
        
    end
end

%legend('0-10 Hz','10-20 Hz','20-30 Hz','30-40 Hz','40-50 Hz');
end