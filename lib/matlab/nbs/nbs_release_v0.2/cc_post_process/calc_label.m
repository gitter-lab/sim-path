function [hcT] = calc_label(c_dmat,knumi,linkF,maskOut)
% Cluster data using hc
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
   
   if (sum(diag(c_dmat(znotempty,znotempty))) > 5) % Some small value
       fprintf(1,'Notify: Transforming into distance metric\n');
       tdist = 1-c_dmat(znotempty,znotempty);
   else 
       tdist = c_dmat(znotempty,znotempty);
   end
   
   if ( any(diag(tdist) ~= 0) )
       fprintf(1,'Warn: nnz on diagnol - Not a distance metric\n');
       tdist = tdist - diag(diag(tdist));
   end
   
   try 
       squareD = squareform(tdist);
   catch e
       fprintf(1,'Error: not clustered due to problem with linkage matrix\n');
       disp(e.message);
       
       hcT = nan(length(c_dmat),1);
       return;
   end
   
   Z = linkage(squareD,linkF);
   hcT = cluster(Z,'maxclust',knumi);

%    if (exist('calcCoph','var') && calcCoph == 1)
%        cophD = cophenet(Z,squareD);
%    end
   
   if (exist('maskOut','var') && maskOut == 1)
       hcT_tmp = hcT;
       hcT = nan(length(c_dmat),1);
       hcT(znotempty) = hcT_tmp;
   end
     
end
