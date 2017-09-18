function [expgeno,cover_fraction,concat_id] = networkPropgateEntrez(netobj,genoEntrez,genoEntrezKey,alpha,isInclude,doNorm,externalKMat)
%
% Network progate a signal on a network
%
%

    zoption.tof = 1e-12;
    zoption.max_iteration = 1000;
    zoption.outlevel = 5;
    
    if (size(genoEntrez,2) ~= length(genoEntrezKey))
        error('Key missmatch error');
    end
            
    [igeno,unMappedGenoMat] = entrezMatToKeyMat(netobj,genoEntrez,genoEntrezKey);
    expgeno = igeno;
    
    cover_fraction = 1 - sum(unMappedGenoMat,2)./sum(genoEntrez,2);
    
    itcnt = 1;             
    if (exist('externalKMat','var') ~= 0 )
        fprintf(1,'External propagation matrix\n');
        adjM = full(externalKMat);
    else        
        adjM = full(netobj.adj_mat_norm);
    end
    
    
    if (exist('doNorm','var') == 1 && doNorm == 1)
        
        expgeno = bsxfun(@times,expgeno,max(sum(expgeno,2),eps).^-1);
  
    end
    
    epgeno_prev = igeno;    
    while (itcnt < zoption.max_iteration)
        
        expgeno = alpha*(mtimesx(expgeno,adjM,'SPEEDOMP')) + (1-alpha)*igeno;
        
       
        l1residual = norm(epgeno_prev - expgeno,1);
        epgeno_prev = expgeno;
        
        if (l1residual < zoption.tof)
            fprintf(1,'Converged %d\n',itcnt);
            
            break;
        end
        itcnt = itcnt + 1;        
    end
    
    if (exist('isInclude','var') == 1 && isInclude == 1)
        gunmapped = sum(unMappedGenoMat) > 1;
        fprintf(1,'Concat %d unmapped gene values\n',sum(gunmapped));        
        % expgeno = [expgeno mean(expgeno(:))*double(genoEntrez(:,gunmapped)>0) ];        
        % zfact = median(max(expgeno));
        zfact = mean(expgeno(:));
        expgeno = [expgeno zfact*double(genoEntrez(:,gunmapped)>0) ];        
        % expgeno = [logit(expgeno,1) double(genoEntrez(:,gunmapped)>0) ];
        concat_id = genoEntrezKey(gunmapped);
    end
    %expgeno = expgeno';
end