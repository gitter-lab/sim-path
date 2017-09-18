function [H,S,numIter,tElapsed,finalResidual,resVal]=nmfrule_symmetric_ding(W,k,option)
% Symmetric NMF with missing values
%
% References:
% Tao Li, Ding C., Jordan MI, Solving Consensus and Semi-supervised Clustering Problems Using Nonnegative Matrix Factorization
%
% This code was extended from Yifeng Li's NMF toolbox 
% Y. Li and A. Ngom, "The non-negative matrix factorization toolbox for 
%      biological data mining," BMC Source Code for Biology and Medicine,
%      vol. 8, pp. 10, April, 2013.
%
% Matan Hofree (mhofre@cs.ucsd.edu)
% V. 0.2.0 (2/15/2013)
%
% This software is Copyright © 2013 
% The Regents of the University of California. All Rights Reserved.
%
% Please see license.txt for full details.



tStart=tic;
optionDefault.iter=1000;
optionDefault.dis=false;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
optionDefault.beta=0.1;
optionDefault.qthrehold = 0.25;

if nargin<3
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

resVal = zeros(optionDefault.iter,1);

SS=isnan(W);
W(SS)=0;
SS=~SS;

betaF = optionDefault.beta;
% iter: number of iterations
[r,c]=size(W); % c is # of samples, r is # of features
H=rand(c,k);
% Y(Y<eps)=0;
H=max(H,eps);
S=H\(W/H');
% A(A<eps)=0;
S=max(S,eps);
XfitPrevious=Inf;
cnt = 1;
for i=1:option.iter

%     HtH = H'*H;
%     S = S.*sqrt((H'*W*H)./(HtH*S*HtH));
%     S = max(S,eps);
        
    HtH = mtimesx(H,'t',H,'SPEEDOMP');
    S = S.*sqrt(mtimesx(H,'t',mtimesx(W,H,'SPEEDOMP'),'SPEEDOMP')./(eps + mtimesx(mtimesx(HtH,S,'SPEEDOMP'),HtH,'SPEEDOMP') ));
    % S = max(S,eps);
    
    
%     WHS = W*H*S;
%     H = H.*sqrt((WHS)./(H*H'*WHS));
%     H = max(H,eps);
    
    WHS = mtimesx(SS.*W,mtimesx(H,S,'SPEEDOMP'),'SPEEDOMP');
    H = H.*sqrt(WHS./mtimesx(mtimesx(H,H,'t','SPEEDOMP').*SS,WHS,'SPEEDOMP'));
    H = max(H,eps);
    
    
    
    if mod(i,10)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        XfitThis=H*S*H';
        fitRes=matrixNorm(XfitPrevious-XfitThis);
        XfitPrevious=XfitThis;
        curRes=norm(W-XfitThis,'fro');
        
        resVal(cnt) = curRes;
        cnt = cnt + 1;
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            if (option.dis)
                s=sprintf('Mutiple update rules based NMF successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,curRes);
                disp(s);
            end
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);

resVal = resVal(1:(cnt-1));
end
