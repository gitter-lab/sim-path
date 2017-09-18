function [hcT,Hout] = calc_label_SNMF(c_dmat,knumi,linkF,maskOut,calcCoph)
% Cluster data using symmetric NMF
% 
%
% Matan Hofree (mhofre@cs.ucsd.edu)
% V. 0.2.0 (2/15/2013)
%
% This software is Copyright © 2013 
% The Regents of the University of California. All Rights Reserved.
%
% Please see license.txt for full details.


   CC_label = [];
   hcT = [];
   if (isempty(c_dmat) || iscell(c_dmat))
       return;
   end
   
   if (exist('maskOut','var') && maskOut == 1)
       znotempty = any(c_dmat ~= 0) & sum(isnan(c_dmat)) < length(c_dmat)/2 & sum(isinf(c_dmat)) < length(c_dmat)/2; 
   else
       znotempty = true(length(c_dmat),1);
       
       if (sum(~any(c_dmat))>0)
           fprintf(1,'Empty rows in cc mat\n');
       end
       
   end
   
   % Fix up similarity matrix
   tdist = c_dmat(znotempty,znotempty);
   tdistR = tdist - diag(diag(tdist)) + eye(size(tdist));
   % [~,tdistR] = sim_to_quantrem(tdistR,0.1);
   % [~,tdistR] = sim_to_midrem(tdistR,0.25);
   % [~,tdistR] = sim_to_knn(tdistR,7);
   
   % [H,S] =  nmfrule_symmetric_dingnan(tdistR.^2,knumi);
   if (any(isnan(tdistR(:))) )       
       [H,S,numIter,tElapsed,finalResidual,resVal] =  nmfrule_symmetric_dingnan(tdistR,knumi);
   else
       [H,S,numIter,tElapsed,finalResidual,resVal] =  nmfrule_symmetric_ding(tdistR,knumi);
   end
   hcT = NMFCluster(H');
   
   
   if (exist('maskOut','var') && maskOut == 1)
       hcT_tmp = hcT;
       hcT = nan(length(c_dmat),1);       
       hcT(znotempty) = hcT_tmp;
       
       Hout = nan(length(c_dmat),size(H,2));
       Hout(znotempty,:) = H;
   end
     
end
