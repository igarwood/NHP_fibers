function fig = plotBeta(beta_a,beta_b)
% Plot frequency and state specific frequency bands
% Beta distribution: p(y) = (y^(a-1)*(1-y)^(b-1))/Beta(a,b)
% Input: for H frequency bands and K states,
%        beta_a: KxH matrix of a 
%        beta_b: KxH matrix of b
     

% Number of states:
K = size(beta_a,1);

% Number of frequency bands:
H = size(beta_a,2);

% Default colors:
defaults = [ [0.6350, 0.0780, 0.1840];...	
            [0.4940, 0.1840, 0.5560];...
            [0.0000, 0.4470, 0.7410];...	          	
          	[0.8500, 0.3250, 0.0980];...	 	          	
          	[0.9290, 0.6940, 0.1250];...	 		 	          	
          	[0.4660, 0.6740, 0.1880];...	 	          	
          	[0.3010, 0.7450, 0.9330]...	 
        ];	
ncolor = size(defaults,1);

defaultLine = max(zeros(size(defaults)),defaults);
for h = 1:H   
    for k = 1:K
        subplot(1,H,h)
        pd = makedist('Beta','a',beta_a(k,h),'b',beta_b(k,h));
        x = linspace(0,1,10000);
        Y = pdf(pd,x);
        plot(x,Y,'Color',defaultLine(mod(k,ncolor)+1,:),'linewidth',2);
        set(gca,'TickLength',[0 0])
        set(gca,'FontSize',16)
        set(gcf,'renderer','Painters')
        hold on
        ax = gca;
        ax.YColor = 'k';
        ax = gca;
        ylabel('f(Y)')
        xlim([0,1])
        box off
        xlabel('Y')        
    end
end


fig = gcf;
end
