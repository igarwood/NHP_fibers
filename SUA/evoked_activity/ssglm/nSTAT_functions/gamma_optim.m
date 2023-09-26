function gradQ = gamma_optim(gamma,basisMat,W_K,x_K,HkAll,dN)
K = length(HkAll);
gradQ=zeros(size(gamma,2),1);
for k=1:K
    Hk=HkAll{k};
    Wk = basisMat*diag(W_K(:,:,k));
    stimK=basisMat*(x_K(:,k));
    histEffect=exp(gamma*Hk')';
    stimEffect = exp(stimK);
    stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
    lambdaDelta = stimEffect.*histEffect;

    gradQ = gradQ + Hk'*dN(k,:)' - Hk'*lambdaDelta;

end
