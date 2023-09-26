function ci = binomial_ci_bs(x)
% x is series of 1s and 0s where sum(x) = binomial rv
% (same function works for any data distribution)

n_bs = 10000;

N = length(x);
p_bs = zeros(1,n_bs);
try
    for n = 1:n_bs
        x_bs = randsample(x,N,1);
        p_bs(n) = sum(x_bs)/N;
    end

    ci(1) = quantile(p_bs,0.025);
    ci(2) = quantile(p_bs,0.975);
catch
    warning('no sample given for binomial_ci_bs');
    ci=zeros(1,2);
end

end