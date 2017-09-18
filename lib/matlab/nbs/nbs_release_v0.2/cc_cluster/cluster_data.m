function [Tnet,H,outlierVect] = cluster_data(gind_sample,gene_id_sample,znetwork,run_param,knnGlap)
       
    Tnet = nan(size(gind_sample,1),length(run_param.K));         
    outlierVect = [];
    
    if ( run_param.isNetwork == 0  || isempty(znetwork) )                
        if (run_param.dis)
            fprintf(1,'No propagation\n');
        end
        expgeno = gind_sample;            
    else        
        if (znetwork.propVal > 0)
            if (run_param.dis)
                fprintf(1,'Using quick propagate\n');
            end
            expgeno = znetwork.networkPropgateEntrezQuick(gind_sample,gene_id_sample,run_param.propV,0,0);
        else
            expgeno = znetwork.networkPropgateEntrez(gind_sample,gene_id_sample,run_param.propV,0,0);
        end
        
        if (run_param.normalize_rows > 0 )
            fprintf(1,'Normalize\n');
            if (exist('intQuantilenorm','file'))
                expgeno = intQuantilenorm(expgeno')';
            else
                expgeno = quantilenorm(expgeno')';
            end
                
        end
    end
          
    if (size(gind_sample,1) < run_param.min_indiv)
        fprintf(1,'Error: Sample too small');
        return;                
    end

    cnt = 0;    
    if (strcmp(run_param.nmf_type,'hc-av-euclid') == 1)
        fprintf(1,'Warn: Using euclid distance');
        Z = linkage(expgeno,'average','euclid');
    elseif (strcmp(run_param.nmf_type,'hc-av-jaccard') == 1)
        fprintf(1,'Using jaccard distance');
        Z = linkage(expgeno,'average','jaccard');
    elseif (strcmp(run_param.nmf_type,'hc-av-zscos') == 1)   
        fprintf(1,'Using zscore cos distance distance');
        Z = linkage(zscore(expgeno,[],2),'average','cos');
    end
    
    for cnum = run_param.K
        cnt = cnt + 1;
        % Network NMF    
        if (strcmp(run_param.nmf_type,'nmf') == 1)
            fprintf(1,'Regular NMF\n');
            [~,H{cnum}] = nmfnnls(expgeno',cnum,run_param.zoptions);            
            indClust = NMFCluster(H{cnum});
        elseif (strcmp(run_param.nmf_type,'nmfOD') == 1)
            fprintf(1,'NMF with outlier detection\n');
            [~,H{cnum},outlierVect{cnum}] = nmfnnls_od(expgeno',cnum,run_param.zoptions);            
            indClust = NMFCluster(H{cnum});            
        elseif (strcmp(run_param.nmf_type,'netnmf') == 1) 
            fprintf(1,'Net NMF\n');
            if (isempty(knnGlap))
                error('Knn glap is empty');
            end
            [~,H{cnum}] =  nbs_nmf_cluster_network_nmf(expgeno',cnum,0,-1,run_param.zoptions.gamma,sparse(knnGlap),run_param.zoptions);
            indClust = NMFCluster(H{cnum});
        elseif (strcmp(run_param.nmf_type(1:2),'hc') == 1)
            fprintf(1,'HC\n');
            indClust = cluster(Z,'maxclust',cnum);
            H{cnum} = clustToAssignmentMatrix(indClust,cnum); 
        else
            error('Cluster option unknown');
        end
        
        Tnet(:,cnt)=indClust;      
    end
    
end
