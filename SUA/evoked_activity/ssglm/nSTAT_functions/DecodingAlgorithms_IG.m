classdef DecodingAlgorithms_IG
% DECODINGALGORITHMS A class that contains static functions for 
% decoding the hidden states of linear discrete stochastic systems or 
% hybrid linear discrete stochastic systems subject to gaussian noise. 
% The observations can come from either a gaussian observation model 
% or via a point process observation model.

% To use, download and install nSTAT toolbox available at:
% https://github.com/iahncajigas/nSTAT
%   
% Modifications to selected functions made by Garwood, I.C., 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original version:
% nSTAT v1 Copyright (C) 2012 Masschusetts Institute of Technology
% Cajigas, I, Malik, WQ, Brown, EN
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% 

    properties
    end
    
    methods (Static)


%% Functions for Point Process State Space Expectation Maximization edited by Indie Garwood
        % PPSS_EMFB implements a forward and backward PPSS_EM. Because the way that the algorithm is setup,
        % we can analyze the data from the first trial to the last (forward) and from
        % the last trial to the first (backward). This approach yields
        % better estimates of the underlying firing rates
        %
        % This version can handle missing trials
        
        function [xKFinal,WKFinal, WkuFinal,Qhat,gammahat,fitResults,stimulus,stimCIs,logll,QhatAll,gammahatAll,nIter]=PPSS_EMFB_missing(A,Q0,x0,dN,missing_trials,fitType,delta,gamma0,windowTimes, numBasis,neuronName)
    %         if(nargin<10 || isempty(neuronName))
    %             neuronName = 1;
    %         end
            dLikelihood(1)=inf;
            if(numel(Q0)==length(Q0)^2)
                Q0=diag(Q0); %turn Q into a vector
            end

            Qhat=Q0;
            gammahat=gamma0;
            xK0=x0;
            cnt=1; tol=1e-2; maxIter=2;%2e3;
            tolAbs = 1e-3;
            tolRel = 1e-2;
            llTol  = 1e-3;
            stoppingCriteria=0;

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    HkAll{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    HkAll{k} = zeros(size(dN,2),1);
                end
                gamma0=0;
            end


            HkAllR=HkAll(end:-1:1);
    %         if(~isempty(windowTimes))
    %             histObj = History(windowTimes,minTime,maxTime);
    %             for k=K:-11:1
    %                 nstr{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    %                 nstr{k}.setMinTime(minTime);
    %                 nstr{k}.setMaxTime(maxTime);
    %                 HkAllR{k} = histObj.computeHistory(nstr{k}).dataToMatrix;
    %             end
    %         else
    %             for k=1:K
    %                 HkAllR{k} = 0;
    %             end
    %             gammahat=0;
    %         end
            Wku = cell(K,K);
            n_backwards_stopped = 0;
            while(stoppingCriteria~=1 && cnt<=maxIter)
                display(['EMFB cycle ' , num2str(cnt),': Forward EM 1']);
                tic
                [xK,WK, Wku,Qhat(:,cnt+1),gammahat(cnt+1,:),logll(cnt),~,~,nIter1,negLL]=DecodingAlgorithms.PPSS_EM_missing(A,Qhat(:,cnt),xK0,dN,missing_trials,fitType,delta,gammahat(cnt,:),windowTimes, numBasis,HkAll,Wku);    
                t = toc;
                display(['  Forward EM 1 completed in ', num2str(t/60),...
                    ' minutes.']);
                if(~negLL) && n_backwards_stopped < 0 % Don't run backward EM if it has failed multiple times in a row
                    display(['EMFB cycle ' , num2str(cnt),': Backward EM']);
                    tic
                    [xKR,~, ~,QhatR(:,cnt+1),gammahatR(cnt+1,:),logllR(cnt),~,~,nIter2,negLL]=DecodingAlgorithms.PPSS_EM_missing(A,Qhat(:,cnt+1),xK(:,end),flipud(dN),flipud(missing_trials),fitType,delta,gammahat(cnt+1,:),windowTimes, numBasis,HkAllR);                    
                    t = toc;
                    display(['  Backward EM completed in ', num2str(t/60),...
                        ' minutes.']);
                    if(~negLL)
                        display(['EMFB cycle ' , num2str(cnt),': Forward EM 2']);
                        tic
                        [xK2,WK2, Wku2,Qhat2,gammahat2,logll2,~,~,nIter3,negLL2]=DecodingAlgorithms.PPSS_EM_missing(A,QhatR(:,cnt+1),xKR(:,end),dN,missing_trials,fitType,delta,gammahatR(cnt+1,:),windowTimes, numBasis,HkAll);
                        t = toc;
                        display(['  Forward EM 2 completed in ',...
                            num2str(t/60), ' minutes.']);
                        if(~negLL2)
                            xK=xK2;
                            WK=WK2;
                            Wku=Wku2;
                           % Wku_indices=Wku_indices2;
                            Qhat(:,cnt+1) = Qhat2;
                            gammahat(cnt+1,:) = gammahat2;
                            logll(cnt) = logll2;
                            n_backwards_stopped = 0;
                        else 
                            n_backwards_stopped = n_backwards_stopped+1;
                        end
                    else
                         n_backwards_stopped = n_backwards_stopped+1;
                    end
                        
                end


                xK0=xK(:,1);
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(logll(cnt)-logll(cnt-1));%./abs(logll(cnt-1));
                end
                cnt=cnt+1;

    %             figure(1)
    %         
    %             subplot(1,2,1); surf(xK);
    %             subplot(1,2,2); plot(logll); ylabel('Log Likelihood');

                dQvals = abs(sqrt(Qhat(:,cnt))-sqrt(Qhat(:,cnt-1)));
                dGamma = abs(gammahat(cnt,:)-gammahat(cnt-1,:));
                dMax = max([dQvals',dGamma]);

                dQRel = max(abs(dQvals./sqrt(Qhat(:,cnt-1))));
                dGammaRel = max(abs(dGamma./gammahat(cnt-1,:)));
                dMaxRel = max([dQRel,dGammaRel]);
    %             dMax
    %             dMaxRel
                if(dMax<tolAbs && dMaxRel<tolRel) & cnt>2
                    stoppingCriteria=1;
                    display(['EMFB converged at iteration:' num2str(cnt) ' b/c change in params was within criteria']);
                end
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0) & cnt>2
                    stoppingCriteria=1;
                    display(['EMFB stopped at iteration:' num2str(cnt) ' b/c change in likelihood was negative']);
                end

            end

            maxLLIndex = find(logll == max(logll),1,'first');
            if(maxLLIndex==1)
                maxLLIndex=cnt-1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
            end

            xKFinal = xK;
            x0Final=xK(:,1);
            WKFinal = WK;
            WkuFinal = Wku;
            QhatAll =Qhat(:,1:maxLLIndex+1);
            Qhat = Qhat(:,maxLLIndex+1);
            gammahatAll =gammahat(1:maxLLIndex+1);
            gammahat = gammahat(maxLLIndex+1,:);
            logll = logll(maxLLIndex);

            K=size(dN,1);
            SumXkTermsFinal = diag(Qhat(:,:,end))*K;
            logllFinal=logll(end);
            McInfo=100;
            McCI = 3000;

            nIter = [];%[nIter1,nIter2,nIter3];
  
            
            K   = size(dN,1); 
            R=size(xK,1);
            logllobs = logll+R*K*log(2*pi)+K/2*log(det(diag(Qhat)))+ 1/2*trace(diag(Qhat)\SumXkTermsFinal);
            InfoMat = [];
            fitResults = [];
            stimCIs = [];
            stimulus = [];
            % InfoMat = DecodingAlgorithms.estimateInfoMat_missing(fitType,dN,missing_trials,HkAll,A,x0Final,xKFinal,WKFinal,WkuFinal,Qhat,gammahat,windowTimes,SumXkTermsFinal,delta,McInfo);
            % fitResults = DecodingAlgorithms.prepareEMResults_missing(fitType,neuronName,dN,HkAll,xKFinal,WKFinal,Qhat,gammahat,windowTimes,delta,InfoMat,logllobs);
%             [stimCIs, stimulus] = DecodingAlgorithms.ComputeStimulusCIs_missing(fitType,xKFinal,WkuFinal,delta,McCI);
    %             

        end
        function [xKFinal,WKFinal, WkuFinal,Qhat,gammahat,logll,QhatAll,gammahatAll,nIter,negLL]=PPSS_EM_missing(A,Q0,x0,dN,missing_trials,fitType,delta,gamma0,windowTimes, numBasis,Hk,Wku_init)
            if(nargin<9 || isempty(numBasis))
                numBasis = 20;
            end
            if(nargin<8 || isempty(windowTimes))
                if(isempty(gamma0))
                    windowTimes =[];
                else
    %                 numWindows =length(gamma0)+1; 
                    windowTimes = 0:delta:(length(gamma0)+1)*delta;
                end
            end
            if(nargin<7)
                gamma0=[];
            end
            if(nargin<6 || isempty(delta))
                delta = .001;
            end
            if(nargin<5)
                fitType = 'poisson';
            end


            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            K=size(dN,1);




    %         tol = 1e-3; %absolute change;
            tolAbs = 1e-2;
            tolRel = 1e-2;
            llTol  = 1e-2;
            cnt=1;

            maxIter = 50;

            if(numel(Q0)==length(Q0)^2)
                Q0=diag(Q0); %turn Q into a vector
            end
            numToKeep=3; %10 Consider memory allocation
            Qhat = zeros(length(Q0),numToKeep);
            Qhat(:,1)=Q0;
            gammahat=zeros(numToKeep,length(gamma0));
            gammahat(1,:)=gamma0;
%             QhatNew=Q0;
%             gammahatNew(1,:)=gamma0;
            cnt=1;
            dLikelihood(1)=inf;
    %         logll(1)=-inf;
            x0hat = x0;
            negLL=0;
         
            %Forward EM
            stoppingCriteria =0;
%             logllNew= -inf;
            storeIndLast = 1;
            
            Wku{1} = Wku_init;
            

            while(stoppingCriteria~=1 && cnt<=maxIter)
                 storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
                 storeIndP1= mod(cnt,numToKeep)+1;
                 storeIndM1= mod(cnt-2,numToKeep)+1;
%                  [xK{storeInd},WK{storeInd},Wku{storeInd},Wku_indices{storeInd},logll(cnt),SumXkTerms,sumPPll]= ...
%                     DecodingAlgorithms.PPSS_EStep_missing(A,Qhat(:,storeInd),x0hat,dN,missing_trials,Hk,fitType,delta,gammahat(storeInd,:),numBasis);
                
                [xK{storeInd},WK{storeInd},Wku{storeInd},logll(cnt),SumXkTerms,sumPPll]= ...
                    DecodingAlgorithms_IG.PPSS_EStep_missing(A,Qhat(:,storeInd),x0hat,dN,missing_trials,Hk,fitType,delta,gammahat(storeInd,:),numBasis,Wku{storeIndLast});

                 [Qhat(:,storeIndP1),gammahat(storeIndP1,:)] = DecodingAlgorithms_IG.PPSS_MStep_missing(dN,missing_trials,Hk,fitType,xK{storeInd},WK{storeInd},gammahat(storeInd,:),delta,SumXkTerms,windowTimes,cnt);

                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(logll(cnt)-logll(cnt-1));%./abs(logll(cnt-1));
                end
% 
              if cnt ==100
                    scrsz = get(0,'ScreenSize');
                    h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
                    subplot(1,2,1); imagesc(xK{storeInd}'); axis xy;
                    subplot(1,2,2); plot(logll); ylabel('Log Likelihood');
                    figure(h)
                elseif (mod(cnt,100)==0)
                    subplot(1,2,1); surf(xK{storeInd});
                    subplot(1,2,2); plot(logll); ylabel('Log Likelihood');
                    figure(h);
                end
                
                dQvals = abs(sqrt(Qhat(:,storeInd))-sqrt(Qhat(:,storeIndM1)));
                dGamma = abs(gammahat(storeInd,:)-gammahat(storeIndM1,:));
                dMax = max([dQvals',dGamma]);

                dQRel = max(abs(dQvals./sqrt(Qhat(:,storeIndM1))));
                dGammaRel = max(abs(dGamma./gammahat(storeIndM1,:)));
                dMaxRel = max([dQRel,dGammaRel]);

             
                cnt=(cnt+1);
                if(dMax<tolAbs && dMaxRel<tolRel) 
                    stoppingCriteria=1;
                    display(['	EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
                    negLL=0;
                end
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)  
                    stoppingCriteria=1;
                    display(['	EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
                    negLL=1;
                end
                
                storeIndLast = storeInd;    
            end


            maxLLIndex  = find(logll == max(logll),1,'first');
            maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
            if(maxLLIndex==1)
%                 maxLLIndex=cnt-1;
                maxLLIndex =1;
                maxLLIndMod = 1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
               maxLLIndMod = 1;
%             else
%                maxLLIndMod = mod(maxLLIndex,numToKeep); 
               
            end
            nIter   = cnt-1;  
%             maxLLIndMod
            xKFinal = xK{maxLLIndMod};
            WKFinal = WK{maxLLIndMod};
            WkuFinal = Wku{maxLLIndMod};
            %Wku_indices = Wku_indices{maxLLIndMod};
            QhatAll =Qhat(:,1:maxLLIndMod);
            Qhat = Qhat(:,maxLLIndMod);
            gammahatAll =gammahat(1:maxLLIndMod);
            gammahat = gammahat(maxLLIndMod,:);
            logll = logll(maxLLIndex);
           
        end
        
        % Subroutines for the PPSS_EM algorithm
        %[x_K,W_K,Wku_flat,Wku_indices,logll,sumXkTerms,sumPPll]...
        function [x_K,W_K,Wku,logll,sumXkTerms,sumPPll]...
                =PPSS_EStep_missing(A,Q,x0,dN,missing_trials,HkAll,fitType,delta,gamma,numBasis,Wku)


             minTime=0;
             maxTime=(size(dN,2)-1)*delta;


             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end
            if(numel(Q)==length(Q))
                Q=diag(Q); %turn Q into a diagonal matrix
            end
            [K,N]   = size(dN); 
            R=size(basisMat,2);

            x_p     = zeros( size(A,2), K );
            x_u     = zeros( size(A,2), K );
            W_p    = zeros( size(A,2),size(A,2), K);
            W_u    = zeros( size(A,2),size(A,2), K );



            for k=1:K

                if(k==1)
                    x_p(:,k)     = A * x0;
                    W_p(:,:,k)   = Q;
                else
                    x_p(:,k)     = A * x_u(:,k-1);
                    W_p(:,:,k)   = A * W_u(:,:,k-1) * A' + Q;
                end

                 sumValVec=zeros(size(W_p,1),1);
                 sumValMat=zeros(size(W_p,2),size(W_p,2));

                if missing_trials(k) == 0 

                    if(strcmp(fitType,'poisson'))
                        Hk=HkAll{k};
                        Wk = basisMat*diag(W_p(:,:,k));
                        stimK=basisMat*x_p(:,k);

                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK);
                        lambdaDelta =stimEffect.*histEffect;
                        GradLogLD =basisMat;
                        JacobianLogLD = zeros(R,R);
                        GradLD = basisMat.*repmat(lambdaDelta,[1 R]);

                        sumValVec = GradLogLD'*dN(k,:)' - diag(GradLD'*basisMat);
                        sumValMat = GradLD'*basisMat;


                    elseif(strcmp(fitType,'binomial'))
                        Hk=HkAll{k};
                        Wk = basisMat*diag(W_p(:,:,k));
                        stimK=basisMat*x_p(:,k);

                        lambdaDelta=exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));  
                        GradLogLD =basisMat.*(repmat(1-lambdaDelta,[1 R]));
                        JacobianLogLD = basisMat.*repmat(lambdaDelta.*(-1+lambdaDelta),[1 R]);
                        GradLD = basisMat.*(repmat(lambdaDelta.*(1-lambdaDelta),[1 R]));
                        JacobianLD = basisMat.*(repmat(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta.^2),[1 R]));

                        sumValVec = GradLogLD'*dN(k,:)' - diag(GradLD'*basisMat);
                        sumValMat = -diag(JacobianLogLD'*dN(k,:)')+ JacobianLD'*basisMat;

                    end  
                end 


                    invW_u             = eye(size(W_p(:,:,k)))/W_p(:,:,k)+ sumValMat;
                    W_u(:,:,k)       = eye(size(invW_u))/invW_u;%+100*diag(eps*rand(size(W_p,1),1));

                 % Maintain Positive Definiteness
                % Make sure eigenvalues are positive
                [vec,val]=eig(W_u(:,:,k) ); val(val<=0)=eps;
                W_u(:,:,k) =vec*val*vec';
                x_u(:,k)  = x_p(:,k)  + W_u(:,:,k)*(sumValVec);

            end

            %[x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u);
            [x_K, W_K,A_k] = DecodingAlgorithms.kalman_smootherFromFiltered_ig(x_p, W_p, x_u, W_u);
%             for k = 1:K
%                 if sum(sum(abs(W_u(:,:,k)-diag(diag(W_u(:,:,k))))))~=0 || ...
%                     sum(sum(abs(W_K(:,:,k)-diag(diag(W_K(:,:,k))))))~=0 || ...
%                     sum(sum(abs(W_p(:,:,k)-diag(diag(W_p(:,:,k))))))~=0 
%                     disp('Off diagonal element nonzero')
%                 end           
%             end
%             Wku_flat=zeros(R,R,((K^2-K)/2+K));
%             Wku_indices=zeros(2,((K^2-K)/2+K));
%             Tk = zeros(R,R,K-1);
%             for k=1:K
%                 index = K*(k-1)+k-k*(k-1)/2;
%                 Wku_flat(:,:,index)=W_K(:,:,k);
%                 Wku_indices(:,index) = [k;k];
%             end
% 
%             for u=K:-1:2
%                 %tic
%                 
%                 for k=(u-1):-1:1
%                     %Tk(:,:,k)=A;
%                     Tk(:,:,k)=A_k(:,:,k);
% %                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
% 
%                      Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
%                      %Dk(:,:,k)=W_u(:,:,k)*A/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
% 
%                      index_a = K*(k-1)+u-k*(k-1)/2;   
%                      index_b = K*(k)+u-k*(k+1)/2; 
%                      %Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
%                      Wku_flat(:,:,index_a)=Dk(:,:,k)*Wku_flat(:,:,index_b);
%                      Wku_indices(:,index_a) = [k;u];
%                 end
%                 %toc
%             end

            % This step requires a large amount of memory; should see if I
            % can use sparse to solve this since the matrix is sym.
           % Wku= cell(K,K);
            
            Tk = zeros(R,R,K-1);
            for k=1:K
                Wku{k,k}=W_K(:,:,k);
            end

% %             The only time you need the full Wku matrix is when
% %             estimating confidence intervals for an estimated model. It
% %             is very slow, so only estimate the terms needed to compute
% %             the LL
%             for u=K:-1:2
%                 for k=(u-1):-1:1
%                     Tk(:,:,k)=A;
%                     %Tk(:,:,k)=A_k(:,:,k);
% %                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
%                      Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
%                     Wku{k,u}=Dk(:,:,k)*Wku{k+1,u};
%                 end
%             end
            
            
           for u=K:-1:2
               k = u-1;
%                 for k=(u-1):-1:1
                    Tk(:,:,k)=A;
                    %Tk(:,:,k)=A_k(:,:,k);
%                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
                    Wku{k,u}=Dk(:,:,k)*Wku{k+1,u};
%                 end
            end
%             Wku=sparse(zeros(R,R,K,K));
%             Tk = zeros(R,R,K-1);
%             for k=1:K
%                 Wku(:,:,k,k)=W_K(:,:,k);
%             end
% 
%             for u=K:-1:2
%                 for k=(u-1):-1:1
%                     Tk(:,:,k)=A;
%                     %Tk(:,:,k)=A_k(:,:,k);
% %                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
%                      Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
%                     Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
%                     Wku(:,:,u,k)=Wku(:,:,k,u);
%                 end
%             end
            
            %All terms
            Sxkxkp1 = zeros(R,R);
            Sxkp1xkp1 = zeros(R,R);
            Sxkxk = zeros(R,R);
            for k=1:K-1
%                Sxkxkp1 = Sxkxkp1+x_u(:,k)*x_K(:,k+1)'+ ...
%                    Lk(:,:,k)*(W_K(:,:,k+1)+(x_K(:,k+1)-x_p(:,k+1))*x_K(:,k+1)');
               u = k+1;
               Sxkxkp1 = Sxkxkp1+Wku{k,k+1}+x_K(:,k)*x_K(:,k+1)';
               %Sxkxkp1 = Sxkxkp1+Wku(:,:,k,k+1)+x_K(:,k)*x_K(:,k+1)';
%                index = K*(k-1)+u-k*(k-1)/2;
%                Sxkxkp1 = Sxkxkp1+Wku_flat(:,:,index)+x_K(:,k)*x_K(:,k+1)';
               Sxkp1xkp1 = Sxkp1xkp1+W_K(:,:,k+1)+x_K(:,k+1)*x_K(:,k+1)';
               Sxkxk = Sxkxk+W_K(:,:,k)+x_K(:,k)*x_K(:,k)';

            end

            sumXkTerms =  Sxkp1xkp1-A*Sxkxkp1-Sxkxkp1'*A'+A*Sxkxk*A'+ ...
                          W_K(:,:,1)+x_K(:,1)*x_K(:,1)' + ... %expected value of xK(1)^2
                          -A*x0*x_K(:,1)' -x_K(:,1)*x0'*A' +A*(x0*x0')*A';


            if(strcmp(fitType,'poisson'))
                sumPPll=0;
                for k=1:K
                    if missing_trials(k) == 0
                        Hk=HkAll{k};
                        Wk = basisMat*diag(W_K(:,:,k));
                        stimK=basisMat*x_K(:,k);
                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
            %             stimEffect=exp(stimK  + Wk*0.5);
                        ExplambdaDelta =stimEffect.*histEffect;
                        ExplogLD = (stimK + (gamma*Hk')');
                        %sumPPll=sum(dN(k,:)'.*ExplogLD - ExplambdaDelta);
                        % IG: above is incorrect; corrected below
                        sumPPll=sumPPll+sum(dN(k,:)'.*ExplogLD - ExplambdaDelta);
                    end
                end
            elseif(strcmp(fitType,'binomial'))

                sumPPll=0;
                for k=1:K
                    Hk=HkAll{k};
                    Wk = basisMat*diag(W_K(:,:,k));
                    stimK=basisMat*x_K(:,k);
                    lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                    ExplambdaDelta=lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;  
            %       logLD = stimK+(gamma*Hk')' - log(1+lambdaDelta)
                    ExplogLD = stimK+(gamma*Hk')' - log(1+exp(stimK+(gamma*Hk')')) -Wk.*(lambdaDelta).*(1-lambdaDelta)*.5;
                    %E(f(x)]=f(x_hat) + 1/2sigma_x^2 * d^2/dx*f(x_hat)
                    %This is applied to log(1+exp(x_K))
            %         
                    sumPPll=sum(dN(k,:)'.*ExplogLD - ExplambdaDelta);
                end

            end
            R=numBasis;
            %logll = -R*K*log(2*pi)-K/2*log(det(Q))  + sumPPll - 1/2*trace(pinv(Q)*sumXkTerms);
            %IG: edited to match publication below
            logll = -R*K/2*log(2*pi)-K/2*log(det(Q))  + sumPPll ...
                - 1/2*trace(pinv(Q)*sumXkTerms);
            
            


        end
        function [Qhat,gamma_new] = PPSS_MStep_missing(dN,missing_trials,HkAll,fitType,x_K,W_K,gamma, delta,sumXkTerms,windowTimes,cnt)
             K=size(dN,1);
             N=size(dN,2);


             sumQ =  diag(diag(sumXkTerms));
             Qhat = sumQ*(1/K);

             [vec,val]=eig(Qhat); val(val<=0)=0.00000001;
             Qhat =vec*val*vec';
             Qhat = (diag(Qhat));
%              if max(Qhat > 0.2)
%                  disp([num2str(max(Qhat)),' ',num2str(cnt)])
%                  Qhat(Qhat>0.2) = 0.2;
%              end

             minTime=0;
             maxTime=(size(dN,2)-1)*delta;

             numBasis = size(x_K,1);
             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end

            options = optimset('Display','off','UseParallel',false,...
                'MaxFunEvals',5000);
%             gamma_new = lsqnonlin(@(gamma) ...
%                 gamma_optim(gamma,basisMat,W_K,x_K,HkAll,dN),...
%                 gamma, [],[],options);
            gamma0 = gamma;
            gamma_new = lsqnonlin(@(gamma) ...
                gamma_optim(gamma,basisMat,W_K(:,:,~missing_trials),...
                x_K(:,~missing_trials),HkAll(~missing_trials),...
                dN(~missing_trials,:)),...
                gamma0, [],[],options);
        end
        function fitResults=prepareEMResults_missing(fitType,neuronNumber,dN,HkAll,xK,WK,Q,gamma,windowTimes,delta,informationMatrix,logll)


            [numBasis, K] =size(xK);
            % Information matrix is very computationally intensive to
            % compute; uncomment this if essential; requires full WKu
            % matrix
            % SE = sqrt(abs(diag(inv(informationMatrix)))); 
            xKbeta = reshape(xK,[numel(xK) 1]);
            seXK=[];
            for k=1:K
                seXK   = [seXK; sqrt(diag(WK(:,:,k)))];
            end
            statsStruct.beta=[xKbeta;(Q(:,end));gamma(end,:)'];
            statsStruct.se  =[seXK];%;SE];
            covarianceLabels = cell(1,numBasis);
            for r=1:numBasis
                if(r<10)
                    covarianceLabels{r} =  ['Q0' num2str(r)];
                else
                    covarianceLabels{r} =  ['Q' num2str(r)];
                end
            end

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end

            nst = cell(1,K);
            if(~isempty(windowTimes))
                histObj{1} = History(windowTimes,minTime,maxTime);
            else
                histObj{1} = [];
            end

            if(isnumeric(neuronNumber))
                name=num2str(neuronNumber);
                if(neuronNumber>0 && neuronNumber<10)
                    name = strcat(num2str(0),name);
                end
                name = ['N' name];  
            else
                name = neuronNumber;
            end

            for k=1:K
                nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta,name);
                nst{k}.setMinTime(minTime);
                nst{k}.setMaxTime(maxTime);

            end

            nCopy = nstColl(nst);
            nCopy = nCopy.toSpikeTrain;
            lambdaData=[];
            cnt=1;

            for k=1:K
                Hk=HkAll{k};
                stimK=basisMat*xK(:,k);


                if(strcmp(fitType,'poisson'))
                    histEffect=exp(gamma(end,:)*Hk')';
                    stimEffect=exp(stimK);
                    lambdaDelta = histEffect.*stimEffect;
                    lambdaData = [lambdaData;lambdaDelta/delta];
                elseif(strcmp(fitType,'binomial'))
                    histEffect=exp(gamma(end,:)*Hk')';
                    stimEffect=exp(stimK);
                    lambdaDelta = histEffect.*stimEffect;
                    lambdaDelta = lambdaDelta./(1+lambdaDelta);
                    lambdaData = [lambdaData;lambdaDelta/delta];
                end


                for r=1:numBasis
                        if(r<10)
                            otherLabels{cnt} = ['b0' num2str(r) '_{' num2str(k) '}']; 
                        else
                            otherLabels{cnt} = ['b' num2str(r) '_{' num2str(k) '}'];
                        end
                        cnt=cnt+1;
                end
            end

            lambdaTime = minTime:delta:(length(lambdaData)-1)*delta;
            nCopy.setMaxTime(max(lambdaTime));
            nCopy.setMinTime(min(lambdaTime));

            numLabels = length(otherLabels);
            if(~isempty(windowTimes))
                histLabels  = histObj{1}.computeHistory(nst{1}).getCovLabelsFromMask;
            else
                histLabels = [];
            end
            otherLabels((numLabels+1):(numLabels+length(covarianceLabels)))=covarianceLabels;
            numLabels = length(otherLabels);

            tc{1} = TrialConfig(otherLabels,sampleRate,histObj,[]); 
            numBasisStr=num2str(numBasis);
            numHistStr = num2str(length(windowTimes)-1);
            if(~isempty(histObj))
                tc{1}.setName(['SSGLM(N_{b}=', numBasisStr,')+Hist(N_{h}=' ,numHistStr,')']);
            else
                tc{1}.setName(['SSGLM(N_{b}=', numBasisStr,')']);
            end
            configColl= ConfigColl(tc);


            otherLabels((numLabels+1):(numLabels+length(histLabels)))=histLabels;




            labels{1}  = otherLabels; % Labels change depending on presence/absense of History or ensCovHist
            if(~isempty(windowTimes))
                numHist{1} = length(histObj{1}.windowTimes)-1;
            else 
                numHist{1}=[];
            end

            ensHistObj{1} = [];
            lambdaIndexStr=1;
            lambda=Covariate(lambdaTime,lambdaData,...
                           '\Lambda(t)','time',...
                           's','Hz',strcat('\lambda_{',lambdaIndexStr,'}'));


            AIC = 2*length(otherLabels)-2*logll;
            BIC = -2*logll+length(otherLabels)*log(length(lambdaData));

            dev=-2*logll;
            b{1} = statsStruct.beta;
            stats{1} = statsStruct;

            distrib{1} =fitType;
            currSpikes=nst;%nspikeColl.getNST(tObj.getNeuronIndFromName(neuronNames));
            for n=1:length(currSpikes)
                currSpikes{n} = currSpikes{n}.nstCopy;
                currSpikes{n}.setName(nCopy.name);
            end
            XvalData{1} = [];
            XvalTime{1} = [];
            spikeTraining = currSpikes;


            fitResults=FitResult(spikeTraining,labels,numHist,histObj,ensHistObj,lambda,b, dev, stats,AIC,BIC,logll,configColl,XvalData,XvalTime,distrib);
            DTCorrection=1;
            makePlot=0;
           Analysis.KSPlot(fitResults,DTCorrection,makePlot);
           Analysis.plotInvGausTrans(fitResults,makePlot);
           Analysis.plotFitResidual(fitResults,[],makePlot); 
        end
        
        function [CIs, stimulus]  = ComputeStimulusCIs_missing(fitType,xK,Wku,delta,Mc,alphaVal)
            if(nargin<7 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<6 ||isempty(Mc))
                Mc=3000;
            end
            [numBasis,K]=size(xK);


           for r=1:numBasis  
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
    %                 stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))/delta;
                    if(strcmp(fitType,'poisson'))
                        stimulusDraw(r,:,c) =  exp(xKDraw(r,:,c))/delta;
                    elseif(strcmp(fitType,'binomial'))
                        stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))./(1+exp(xKDraw(r,:,c)))/delta;
                    end
                end
           end

           CIs = zeros(size(xK,1),size(xK,2),2);
           for r=1:numBasis
               for k=1:K
                   [f,x] = ecdf(squeeze(stimulusDraw(r,k,:)));
                    CIs(r,k,1) = x(find(f<alphaVal/2,1,'last'));
                    CIs(r,k,2) = x(find(f>(1-(alphaVal/2)),1,'first'));
               end
           end

           if(nargout==2)
               if(strcmp(fitType,'poisson'))
                    stimulus =  exp(xK)/delta;
               elseif(strcmp(fitType,'binomial'))
                    stimulus = exp(xK)./(1+exp(xK))/delta;
               end
           end


        end
       
        function [CIs, stimulus]  = ComputeStimulusCIs_group_missing(fitType,xK,Wku,delta,groupK,Mc,alphaVal)
            % Compute stimulus and CIs for equally sized groups of trials
            % (e.g. 10 at the beginning of the session and 10 at the end)
            if(nargin<6 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<5 ||isempty(Mc))
                Mc=3000;
            end
            [numBasis,K]=size(xK);

            numGroups = size(groupK,1); 
            
           for r=1:numBasis  
                
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
    %                 stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))/delta;
                    if(strcmp(fitType,'poisson'))
                        stimulusDraw(r,:,c) =  exp(xKDraw(r,:,c))/delta;
                    elseif(strcmp(fitType,'binomial'))
                        stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))./(1+exp(xKDraw(r,:,c)))/delta;
                    end
                end
           end

           CIs = zeros(size(xK,1),numGroups,2);
           stimulus = zeros(size(xK,1),numGroups,1);
           
           for r=1:numBasis
               for g = 1:numGroups
                   [f,x] = ecdf(reshape(stimulusDraw(r,groupK(g,:),:),size(groupK,2)*Mc,1));
                   CIs(r,g,1) = x(find(f<alphaVal/2,1,'last'));
                   CIs(r,g,2) = x(find(f>(1-(alphaVal/2)),1,'first'));
                   stimulus(r) = x(find(f<0.5,1,'last'));
               end
           end

%            if(nargout==2)
%                if(strcmp(fitType,'poisson'))
%                     stimulus =  exp(xK)/delta;
%                elseif(strcmp(fitType,'binomial'))
%                     stimulus = exp(xK)./(1+exp(xK))/delta;
%                end
%            end


        end
        
        function InfoMatrix=estimateInfoMat_missing(fitType,dN,...
                missing_trials,HkAll,A,x0,xK,WK,Wku,Q,gamma,...
                windowTimes,SumXkTerms,delta,Mc)
            if(nargin<14)
                Mc=500;
            end

            [K,N]=size(dN);
            if(~isempty(windowTimes))
                J=max(size(gamma(end,:)));
            else
                J=0;
            end

            R=size(Q,1);
            numBasis = R;

            % The complete data information matrix
            Ic=zeros(J+R,J+R);
            Q=(diag(Q)); % Make sure Q is diagonal matrix


            X=((SumXkTerms));
            Ic(1:R,1:R) = K/2*eye(size(Q))/Q^2 +X'/Q^3;


            % Compute information of history terms
            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
    %         nst = cell(1,K);
    %         if(~isempty(windowTimes))
    %             histObj = History(windowTimes,minTime,maxTime);
    %             for k=1:K
    %                 nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    %                 nst{k}.setMinTime(minTime);
    %                 nst{k}.setMaxTime(maxTime);
    %                 Hn{k} = histObj.computeHistory(nst{k}).dataToMatrix;
    %             end
    %         else
    %             for k=1:K
    %                 Hn{k} = 0;
    %             end
    %             gamma=0;
    %         end

             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end

            jacQ =zeros(size(gamma,2),size(gamma,2));
            if(strcmp(fitType,'poisson'))
                
                for k=1:K
                    if missing_trials(k) == 0
                        Hk=HkAll{k};

                        Wk = basisMat*diag(WK(:,:,k));
                        stimK=basisMat*(xK(:,k));
                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
                        lambdaDelta = stimEffect.*histEffect;

                        jacQ  = jacQ  - (Hk.*repmat(lambdaDelta,[1 size(Hk,2)]))'*Hk;
                    end
                 end

             elseif(strcmp(fitType,'binomial'))
                 for k=1:K
                    Hk=HkAll{k};
                    Wk = basisMat*diag(WK(:,:,k));
                    stimK=basisMat*(xK(:,k));

                    histEffect=exp(gamma*Hk')';
                    stimEffect=exp(stimK);
                    C = stimEffect.*histEffect;
                    M = 1./C;
                    lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                    ExpLambdaDelta = lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;
                    ExpLDSquaredTimesInvExp = (lambdaDelta).^2.*1./C;
                    ExpLDCubedTimesInvExpSquared = (lambdaDelta).^3.*M.^2 +Wk/2.*(3.*M.^4.*lambdaDelta.^3+12.*lambdaDelta.^3.*M.^3-12.*M.^4.*lambdaDelta.^4);

                    jacQ  = jacQ  - (Hk.*repmat(ExpLDSquaredTimesInvExp.*dN(k,:)',[1,size(Hk,2)]))'*Hk ...
                                  - (Hk.*repmat(ExpLDSquaredTimesInvExp,[1,size(Hk,2)]))'*Hk ...
                                  - (Hk.*repmat(2*ExpLDCubedTimesInvExpSquared,[1,size(Hk,2)]))'*Hk;

                 end


            end           

            Ic(1:R,1:R)=K*eye(size(Q))/(2*(Q)^2)+(eye(size(Q))/((Q)^3))*SumXkTerms;

            if(~isempty(windowTimes))
                Ic((R+1):(R+J),(R+1):(R+J)) = -jacQ;
            end
            xKDraw = zeros(numBasis,K,Mc);
            for r=1:numBasis 
                WkuTemp = zeros(K,K);
                for j = 1:K
                    for k = j:K
                        WkuTemp(j,k) = Wku{j,k}(r,r);
                        WkuTemp(k,j) = Wku{j,k}(r,r);
                    end
                end
                

%                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
            end



            Im=zeros(J+R,J+R);
            ImMC=zeros(J+R,J+R);

            for c=1:Mc

                gradQGammahat=zeros(size(gamma,2),1);
                gradQQhat=zeros(1,R);        
                if(strcmp(fitType,'poisson'))
                    for k=1:K
                        if missing_trials(k) == 0
                            Hk=HkAll{k};
                            stimK=basisMat*(xKDraw(:,k,c));
                            histEffect=exp(gamma*Hk')';
                            stimEffect=exp(stimK);
                            lambdaDelta = stimEffect.*histEffect;
                            gradQGammahat = gradQGammahat + Hk'*dN(k,:)' - Hk'*lambdaDelta;
                            if(k==1)
                                gradQQhat = ((xKDraw(:,k,c)-A*x0).*(xKDraw(:,k,c)-A*x0));
                            else
                                gradQQhat = gradQQhat+((xKDraw(:,k,c)-A*xKDraw(:,k-1,c)).*(xKDraw(:,k,c)-A*xKDraw(:,k-1,c)));
                            end
                        end
                    end
                elseif(strcmp(fitType,'binomial'))
                     for k=1:K
                        Hk=HkAll{k};
                        Wk = basisMat*diag(WK(:,:,k));
                        stimK=basisMat*(xKDraw(:,k,c));

                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK);
    %                   
                        C = stimEffect.*histEffect;
                        M = 1./C;
                        lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                        ExpLambdaDelta = lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;
                        ExpLDSquaredTimesInvExp = (lambdaDelta).^2.*1./C;
                        ExpLDCubedTimesInvExpSquared = (lambdaDelta).^3.*M.^2 +Wk/2.*(3.*M.^4.*lambdaDelta.^3+12.*lambdaDelta.^3.*M.^3-12.*M.^4.*lambdaDelta.^4);


                        gradQGammahat = gradQGammahat + (Hk.*repmat(1-ExpLambdaDelta,[1,size(Hk,2)]))'*dN(k,:)' ...
                                          - (Hk.*repmat(ExpLDSquaredTimesInvExp./lambdaDelta,[1,size(Hk,2)]))'*lambdaDelta;
                        if(k==1)
                            gradQQhat = ((xKDraw(:,k,c)-A*x0).*(xKDraw(:,k,c)-A*x0));
                        else
                            gradQQhat = gradQQhat+((xKDraw(:,k,c)-A*xKDraw(:,k-1,c)).*(xKDraw(:,k,c)-A*xKDraw(:,k-1,c)));
                        end
                     end


                end

                gradQQhat = .5*eye(size(Q))/Q*gradQQhat - diag(K/2*eye(size(Q))/Q^2);
                ImMC(1:R,1:R)=ImMC(1:R,1:R)+gradQQhat*gradQQhat';
                if(~isempty(windowTimes))
                    ImMC((R+1):(R+J),(R+1):(R+J)) = ImMC((R+1):(R+J),(R+1):(R+J))+diag(diag(gradQGammahat*gradQGammahat'));
                end
            end
            Im=ImMC/Mc;

            InfoMatrix=Ic-Im; % Observed information matrix



        end
        function [spikeRateSig, ProbMat,sigMat]=computeSpikeRateCIs_missing(xK,Wku,dN,t0,tf,fitType,delta,gamma,windowTimes,Mc,alphaVal)
             if(nargin<11 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<10 ||isempty(Mc))
                Mc=500;
            end

            [numBasis,K]=size(xK);

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end


    %         K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    Hk{k} = 0;
                end
                gamma=0;
            end

           for r=1:numBasis  
               WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
           end

           time=minTime:delta:maxTime;
           for c=1:Mc
               for k=1:K

                   if(strcmp(fitType,'poisson'))
                        stimK=basisMat*xKDraw(:,k,c);
                        histEffect=exp(gamma*Hk{k}')';
                        stimEffect=exp(stimK);
                        lambdaDelta(:,k,c) =stimEffect.*histEffect;
                   elseif(strcmp(fitType,'binomial'))
                        stimK=basisMat*xKDraw(:,k,c);
                        lambdaDelta(:,k,c)=exp(stimK+(gamma*Hk{k}')')./(1+exp(stimK+(gamma*Hk{k}')'));  
                   end  


               end
               lambdaC=Covariate(time,lambdaDelta(:,:,c)/delta,'\Lambda(t)');
               lambdaCInt= lambdaC.integral;
               spikeRate(c,:) = (1/(tf-t0))*(lambdaCInt.getValueAt(tf)-lambdaCInt.getValueAt(t0));

           end

           CIs = zeros(K,2);
           for k=1:K
               [f,x] = ecdf(spikeRate(:,k));
                CIs(k,1) = x(find(f<alphaVal,1,'last'));
                CIs(k,2) = x(find(f>(1-(alphaVal)),1,'first'));
           end
           spikeRateSig = Covariate(1:K, mean(spikeRate),['(' num2str(tf) '-' num2str(t0) ')^-1 * \Lambda(' num2str(tf) '-' num2str(t0) ')'],'Trial','k','Hz');
           ciSpikeRate = ConfidenceInterval(1:K,CIs,'CI_{spikeRate}','Trial','k','Hz');
           spikeRateSig.setConfInterval(ciSpikeRate);


           if(nargout>1)
               ProbMat = zeros(K,K);
               for k=1:K
                   for m=(k+1):K

                       ProbMat(k,m)=sum(spikeRate(:,m)>spikeRate(:,k))./Mc;
                   end
               end
           end


           if(nargout>2)
                sigMat= double(ProbMat>(1-alphaVal));
           end


        end
        function [spikeRateSig, ProbMat,sigMat]=computeSpikeRateDiffCIs_missing(xK,Wku,dN,time1,time2,fitType,delta,gamma,windowTimes,Mc,alphaVal)
             if(nargin<11 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<10 ||isempty(Mc))
                Mc=500;
            end

            [numBasis,K]=size(xK);

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end


    %         K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    Hk{k} = 0;
                end
                gamma=0;
            end

           for r=1:numBasis  
               WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
           end

           timeWindow=minTime:delta:maxTime;
           for c=1:Mc
               for k=1:K

                   if(strcmp(fitType,'poisson'))
                        stimK=basisMat*xKDraw(:,k,c);
                        histEffect=exp(gamma*Hk{k}')';
                        stimEffect=exp(stimK);
                        lambdaDelta(:,k,c) =stimEffect.*histEffect;
                   elseif(strcmp(fitType,'binomial'))
                        stimK=basisMat*xKDraw(:,k,c);
                        lambdaDelta(:,k,c)=exp(stimK+(gamma*Hk{k}')')./(1+exp(stimK+(gamma*Hk{k}')'));  
                   end  


               end
               lambdaC=Covariate(timeWindow,lambdaDelta(:,:,c)/delta,'\Lambda(t)');
               lambdaCInt= lambdaC.integral;
               spikeRate(c,:) = (1/(max(time1)-min(time1)))*(lambdaCInt.getValueAt(max(time1))-lambdaCInt.getValueAt(min(time1))) ...
                                - (1/(max(time2)-min(time2)))*(lambdaCInt.getValueAt(max(time2))-lambdaCInt.getValueAt(min(time2)));


           end

           CIs = zeros(K,2);
           for k=1:K
               [f,x] = ecdf(spikeRate(:,k));
                CIs(k,1) = x(find(f<alphaVal,1,'last')); %not alpha/2 since this is a once sided comparison
                CIs(k,2) = x(find(f>(1-(alphaVal)),1,'first'));
           end
           spikeRateSig = Covariate(1:K, mean(spikeRate),['(t_{1f}-t_{1o})^-1 * \Lambda(t_{1f}-t_{1o}) - (t_{2f}-t_{2o})^-1 * \Lambda(t_{2f}-t_{2o}) '],'Trial','k','Hz');
           ciSpikeRate = ConfidenceInterval(1:K,CIs,'CI_{spikeRate}','Trial','k','Hz');
           spikeRateSig.setConfInterval(ciSpikeRate);


           if(nargout>1)
               ProbMat = zeros(K,K);
               for k=1:K
                   for m=(k+1):K

                       ProbMat(k,m)=sum(spikeRate(:,m)>spikeRate(:,k))./Mc;
                   end
               end
           end


           if(nargout>2)
                sigMat= double(ProbMat>(1-alphaVal));
           end


        end
    end
end