function [med_ci] = binomial_compare_bs(x,y)
% x and y are series of 1s and 0s where sum(x) = binomial rv
% ci is the difference between probability of x and probability of y.
% (p_x-p_y)
% (same function works for any data distribution)

n_bs = 10000;

N_x = length(x);
N_y = length(y);
p_diff_bs = zeros(1,n_bs);
try
    for n = 1:n_bs
        x_bs = randsample(x,N_x,1);
        y_bs = randsample(y,N_y,1);
        p_x = sum(x_bs)/N_x;
        p_y = sum(y_bs)/N_y;
        p_diff_bs(n) = p_x-p_y;
    end

    med = quantile(p_diff_bs,0.5);
    ci(1) = quantile(p_diff_bs,0.025);
    ci(2) = quantile(p_diff_bs,0.975);

    med_ci = [med,ci];
catch
    warning('empty sample input to binomial_compare_bs');
    med_ci = zeros(1,3);
end

end