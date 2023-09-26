%% Extract behavior statistics for all sessions and plot summary statistics

extract_behavior_stats;

%%
% Default colors for plots:
defaults = [0.6350    0.0780    0.1840;...
        0.0000    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330];  
    
n_session = length(sessions);

figure
y_plot = zeros(1,n_session);
ci_neg = zeros(1,n_session);
ci_pos = zeros(1,n_session);


x_plot = 1:n_session;
for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.accuracy.diff.inhib(1);
    ci_neg(n) = -1*session_stats{n}.accuracy.diff.inhib(2);
    ci_pos(n) = -1*session_stats{n}.accuracy.diff.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

plot([0,n_session+1],[0,0],'k--')

x_plot = x_plot+0.25
for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.accuracy.diff.all.inhib(1);
    ci_neg(n) = -1*session_stats{n}.accuracy.diff.all.inhib(2);
    ci_pos(n) = -1*session_stats{n}.accuracy.diff.all.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end
ylabel('Change in accuracy')
%%
figure
y_plot = zeros(1,n_session);
ci_neg = zeros(1,n_session);
ci_pos = zeros(1,n_session);

for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.complete.diff.inhib(1);
    ci_neg(n) = -1*session_stats{n}.complete.diff.inhib(2);
    ci_pos(n) = -1*session_stats{n}.complete.diff.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

plot([0,n_session+1],[0,0],'k--')

x_plot = x_plot+0.25
for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.complete.diff.all.inhib(1);
    ci_neg(n) = -1*session_stats{n}.complete.diff.all.inhib(2);
    ci_pos(n) = -1*session_stats{n}.complete.diff.all.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

ylabel('Change in trial completion')
%%
figure
y_plot = zeros(1,n_session);
ci_neg = zeros(1,n_session);
ci_pos = zeros(1,n_session);

for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.lr_bias.diff.inhib(1);
    ci_neg(n) = -1*session_stats{n}.lr_bias.diff.inhib(2);
    ci_pos(n) = -1*session_stats{n}.lr_bias.diff.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

plot([0,n_session+1],[0,0],'k--')

x_plot = x_plot+0.25
for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.lr_bias.diff.all.inhib(1);
    ci_neg(n) = -1*session_stats{n}.lr_bias.diff.all.inhib(2);
    ci_pos(n) = -1*session_stats{n}.lr_bias.diff.all.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end
ylabel('Change in left/right bias')



%%
figure
y_plot = zeros(1,n_session);
ci_neg = zeros(1,n_session);
ci_pos = zeros(1,n_session);

for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.reaction_time.diff.inhib(1);
    ci_neg(n) = -1*session_stats{n}.reaction_time.diff.inhib(2);
    ci_pos(n) = -1*session_stats{n}.reaction_time.diff.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

plot([0,n_session+1],[0,0],'k--')

x_plot = x_plot+0.25
for n = 1:n_session
    y_plot(n) = -1*session_stats{n}.reaction_time.diff.all.inhib(1);
    ci_neg(n) = -1*session_stats{n}.reaction_time.diff.all.inhib(2);
    ci_pos(n) = -1*session_stats{n}.reaction_time.diff.all.inhib(3);
    plotConfInterval(x_plot(n),y_plot(n),ci_neg(n), ci_pos(n),...
        defaults(session_tag(n),:));
    hold on
end

ylabel('Change in reaction time (s)')