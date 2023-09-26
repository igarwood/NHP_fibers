function [est_lamb,est_lamb_ci] = estimate_lambda(X,B,stats,n_Y)

n_trials = size(X,2);

est_lamb = zeros(n_Y,n_trials);
est_ci = zeros(n_Y,n_trials,2);

for n = 1:n_trials
    [est_lamb(:,n),dylo,dyhi] = ...
        glmval(B,squeeze(X(:,n,:)),'log',stats,'constant','off');
    est_lamb_ci(:,n,1) = est_lamb(:,n)-dylo;
    est_lamb_ci(:,n,2) = est_lamb(:,n)+dyhi;   
end