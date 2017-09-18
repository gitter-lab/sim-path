function [W,H,tstats] = ...
        nbs_nmf_cluster_network_nmf(Y,k,falpha,fbeta,fgamma,K,option,Hext,Wext)
% NMF based on multiple update rules: Y=WH, s.t. Y,W,H>=0.
% With regularization: 
%           alpha|W| + beta|H| + gamma|WL|
%           Where L is cholesky transform on some kenel matrix (e.g. Graph
%           laplacian)
% 
% Input:
% 	Y - input matrix to factorize
%   k - number of dimensions in factorization
%   falpha - (deprecated)
%   fbeta - (deprecated)
%   fgamma - Initial gamma value 
%   K - A graph Laplacian matrix 
%   options - Runtime option
%   Hext - Initializations for the H matrix
%   Wext - Initializations for the W matrix
% Output:
%   W - A set of cluster prototype vectors
%   H - A set of loadings for the cluster prototypes
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
% optionDefault.distance='ls';
optionDefault.iter=200;
optionDefault.dis=true;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
optionDefault.iterPerformanceSample = 50;
optionDefault.isOffset = 0;
optionDefault.optGamma = 1;
optionDefault.optGammaIterMax = optionDefault.iter/2;
optionDefault.optGammaIterMin = 0;

if exist('option','var') ~= 0 
    option=mergeOption(option,optionDefault);    
else
    option=optionDefault;
end
if ( exist('niter','var') ~= 0 && ~isempty(niter))
    option.iter = niter;
end

resVal = zeros(option.iter,1);
resVal_Kreg = zeros(option.iter,1);
fitResVect = zeros(option.iter,1);
fitGamma = zeros(option.iter,1);

% iter: number of iterations
[r,c]=size(Y); % c is # of samples, r is # of features
if (exist('Hext','var') ~= 1 || isempty(Hext))
    H=rand(k,c);
else
    H = Hext;
    if option.dis
        fprintf(1,'Using prespecified H\n');
    end
end

if (exist('Wext','var') ~= 1 || isempty(Wext))
    H=max(H,eps);
    W=Y/H;
    W = W./repmat(sum(W),length(W),1);   
    W=max(W,eps);
elseif (exist('Hext','var') ~= 1 || isempty(Hext))
    fprintf(1,'Using prespecified W\n');    
    W = Wext;
    H = W\Y;
else     
    fprintf(1,'Using prespecified W and H\n');    
    W = Wext;
    H = Hext;
end

% if ( exist('fgamma_min','var') ~= 0 )
%     min_gamma = fgamma_min;
% else
%     min_gamma = 0;
% end

% fix K size
if ( r > length(K) )
    fprintf(1,'Warn: Added diagnol elements to fix size of KNN matrix\n');
    Ktag = speye(r);
    Ktag(1:length(K),1:length(K)) = K;
    if (issparse(K))
        Ktag = speye(r);
        Ktag(1:length(K),1:length(K)) = K;
        K = sparse(Ktag);
    else
        Ktag = eye(r);
        Ktag(1:length(K),1:length(K)) = K;

        K = Ktag;
    end
end
%% Assumes that the input K is a simple graph laplacian matrix
assert(all(sum(K) < 10e-10),'Error seems input K is not a simple graph laplacian matrix');
Dm = diag(diag(K));
Km = Dm - K;

Ydf = Y;

XfitPrevious=inf;
cnt=1;

for i=1:option.iter         
    
    KWmat_D = mtimesx(Dm,W,'SPEEDOMP');
    KWmat_W = mtimesx(Km,W,'SPEEDOMP');
        
    % Some performance outputs
    if mod(i-1,option.iterPerformanceSample)==0 || i==option.iter
        % Figure out a gamma value
        KWmat = mtimesx(K,W,'SPEEDOMP');
        
        Kres = sqrt(trace(mtimesx(W,'t',KWmat,'SPEEDOMP')));
        XfitThis=mtimesx(W,H,'SPEEDOMP');
        WHres = norm(Ydf - XfitThis,'fro');
        
        if (~isinf(XfitPrevious))
            fitRes=norm(XfitPrevious-XfitThis,'fro');
        else 
            fitRes=nan;
        end
        XfitPrevious=XfitThis;
        
        curRes=WHres;
                
        resVal(cnt) = curRes;        
        resVal_Kreg(cnt) = Kres;
        fitResVect(cnt) = fitRes;
        
        if option.dis
            fprintf(1,'Iterating >>>>>> %dth\tMat-res:%f\tK-res:%f\tSum:%f\tGamma:%f\tWfrob:%f\n',i-1,curRes,resVal_Kreg(cnt),curRes + resVal_Kreg(cnt),fgamma,norm(W,'fro'));
        end
        
        cnt = cnt + 1;
        if (option.tof>=fitRes) || option.residual>=curRes || i==option.iter
            if option.dis
                s=sprintf('Mutiple update rules based NMF successes! \n # of iterations is %0.0d. \n The final residual is %f',i,curRes);
                disp(s);
            end
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
    
    if ((option.optGamma ~= 0) && (i < option.optGammaIterMax) && (i > option.optGammaIterMin))
        fgamma_factor = round((WHres/Kres)*option.optGamma);                        
        
        fgamma = fgamma_factor;
    end
    
    fitGamma(i) = fgamma;
         
    W = W.*((mtimesx(Ydf,H,'t','SPEEDOMP')  + fgamma*KWmat_W + eps )./( mtimesx(W,mtimesx(H,H,'t','SPEEDOMP'),'SPEEDOMP') + fgamma*KWmat_D + eps));    
    
    W = max(W,eps);
    W = bsxfun(@times,W,max(sum(W),eps).^-1);
    
    if (fbeta == 0)
        H=H.*(mtimesx(W,'t',Ydf,'SPEEDOMP')./mtimesx(mtimesx(W,'t',W,'SPEEDOMP'),H,'SPEEDOMP'));
    elseif (fbeta == -1)    
        H=fcnnls(W,Ydf);
%     elseif (falpha ~= 0)
%         H=H.*((mtimesx(W,'t',Ydf,'SPEEDOMP') + falpha*KHmat_W + eps)./(mtimesx(mtimesx(W,'t',W,'SPEEDOMP'),H,'SPEEDOMP') + falpha*KHmat_D + eps));    
%     else 
%         We=[W;sqrt(fbeta)*ones(1,k)];
%         Ydf0=[Ydf;zeros(1,c)];
%         H=fcnnls(We,Ydf0);
    end
    H=max(H,eps);        
end
tElapsed=toc(tStart);
resVal = resVal(1:(cnt-1));
resVal_Kreg = resVal_Kreg(1:(cnt-1));
fitResVect = fitResVect(2:cnt-1);
fitGamma = fitGamma(1:i);

tstats.numIter = numIter;
tstats.tElapsed = tElapsed;
tstats.finalResidual = finalResidual;
tstats.resVal = resVal;
tstats.resVal_Kreg = resVal_Kreg;
tstats.fitResVect = fitResVect;
tstats.fitGamma = fitGamma;

end
