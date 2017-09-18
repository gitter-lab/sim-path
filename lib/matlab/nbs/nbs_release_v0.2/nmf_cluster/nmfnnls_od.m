function [W,H,finalIsOutlier,numIter,tElapsed,finalResidual]=nmfnnls_od(X,k,option)
% NMF with outliers detection and removal
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
%%%%
    
    tStart=tic;
    optionDefault.iter=1000;
    optionDefault.calcResid = 20;
    optionDefault.dis=true;
    optionDefault.residual=1e-4;
    optionDefault.tof=1e-4;
    optionDefault.odquantile = 0.1;
    optionDefault.odminbase = 0.03;
    optionDefault.odCalc = 23;
    optionDefault.odmaxIter = floor(optionDefault.iter*0.33);
    
    if nargin<3
        option=optionDefault;
    else
        option=mergeOption(option,optionDefault);
    end
    
    % iter: number of iterations
    [r,c]=size(X); % c is # of samples, r is # of features
    H=rand(k,c);
    XfitPrevious=Inf;
    
    odnumberSample = floor(c*option.odquantile);
    odminBaseSupport = floor(c*option.odminbase);
    Xtag = X;
    isOutlierBaseFixed = 1;
    
    for i=1:option.iter
        W=kfcnnls(H.',Xtag.');
        W=W.';
        W=normc(W);
        H=kfcnnls(W,Xtag);
        
        nonOutlierRound = 1;
        
        % Outlier detection
        % Do this every odCalc upto odmaxIter and one last time just qunatile
        % outliers
        if mod(i,option.odCalc)==0 && (i < option.odmaxIter || isOutlierBaseFixed ~= 0)
            nonOutlierRound = 0;
            
            % Detect and remov a quantile of outlier samples
            if ( isOutlierBaseFixed ~= 0 )
                isOutlierBaseFixed = 0; % This is a quantile fix round
                
                Hfull=kfcnnls(W,X);
                XfitPrevious = W*Hfull;
                
                XfitResid = sum((X - XfitPrevious).^2);
                [~,idx] = sort(XfitResid,'Descend');
                
                if (option.dis)
                    fprintf(1,'Dropping outliers:');
                    fprintf(1,' %d',idx(1:odnumberSample));
                    fprintf(1,'\n');
                end
                
                nonOutlierSet = true(c,1);
                nonOutlierSet(idx(1:odnumberSample)) = 0;
                
                Xtag = X(:,nonOutlierSet);
                H = Hfull(:,nonOutlierSet);
                XfitPrevious = XfitPrevious(:,nonOutlierSet);
            else %  Detect and remove outlier base vectors
                isOutlierBaseFixed = 1;
                
                [~,idx] = max(H);
                Hassign = clustToAssignmentMatrix(idx,k);
                countMax = sum(Hassign,2)';
                
                if any(countMax < odminBaseSupport)
                    [~,outlierBaseIdx] = min(countMax);
                    [~,outlierMaxSplit] = max(countMax);
                    
                    odBaseDrop = any(Hassign(outlierBaseIdx,:),1);
                    Xtag = Xtag(:,~odBaseDrop);
                    H = H(:,~odBaseDrop);
                    
                    % H(outlierBaseIdx,:) = rand(nnz(outlierBaseIdx),size(H,2));
                    % Try and find a good split into the removed vector
                    %
                    H = split_H_matrix(H,outlierBaseIdx,outlierMaxSplit);
                    
                    XfitPrevious=XfitPrevious(:,~odBaseDrop);
                    
                    if (option.dis)
                        fprintf(1,'Dropped base vectors idx %d - (%d samples)\n',outlierBaseIdx,sum(odBaseDrop));
                    end
                end
            end
        end
        
        
        if mod(i,option.calcResid)==0 || i==option.iter
            if option.dis
                disp(['Iterating >>>>>> ', num2str(i),'th']);
            end
            XfitThis=W*H;
            fitRes=matrixNorm(XfitPrevious-XfitThis);
            XfitPrevious=XfitThis;
            curRes=norm(Xtag-XfitThis,'fro');
            if (option.tof>=fitRes && nonOutlierRound) || option.residual>=curRes || i==option.iter
                if option.dis
                    s=sprintf('NNLS based NMF successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,curRes);
                    disp(s);
                end
                numIter=i;
                finalResidual=curRes;
                break;
            end
        end
    end
    % Fix final H to include outliers
    H=kfcnnls(W,X);
    
    XfitPrevious = W*H;
    
    XfitResid = sum((X - XfitPrevious).^2);
    [~,idx] = sort(XfitResid,'Descend');
    
    finalIsOutlier = true(c,1);
    finalIsOutlier(idx(1:odnumberSample)) = 0;
    
    tElapsed=toc(tStart);
end

function H = split_H_matrix(H,outlierBaseIdx,outlierMaxSplit)
    %
    % Cluster the loaddings to derive a new split of the bigest group
    %
    ccoptions.dis = 0;
    
    [~,idx] = max(H);
    Hassign = clustToAssignmentMatrix(idx,size(H,1));
    
    Hselect = logical(Hassign(outlierMaxSplit,:));
    
    [~,zsplit] = nmfnnls(H(:,Hselect),2,ccoptions);
    [~,idx] = max(zsplit);
    
    HselectPos = find(Hselect);
    
    H(outlierBaseIdx,:) = 0; % rand(1,length(H));
    H(outlierBaseIdx,HselectPos(idx==1)) = H(outlierMaxSplit,HselectPos(idx==1));
    H(outlierMaxSplit,HselectPos(idx==1)) = 0; % rand(1,sum(idx==1));
end
