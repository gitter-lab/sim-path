function [cc_clust] = hmat_to_cc_cluster(cc_type,kvect,hard_label_vect,H,Nsample_vect,Nind,outlierVect)
    
   nsample = length(H);
   kvect = kvect(:)';

   if (exist('outlierVect','var') == 0)
       outlierVect = [];
   end
   
   % cc_clust.clust_hard = [];
   % cc_clust.clust_soft = [];
   
   % If we do not have the hard vect results we extract these from the H 
   if (isempty(hard_label_vect) && any(strcmp(cc_type,'hard')) )
       hard_label_vect = cell(nsample,1);
       for i = 1:nsample
           
           hard_label_vect{i} = nan(size(H{i}{end},2),length(kvect));
           cnt = 1;
           for cnum = kvect
               hard_label_vect{i}(:,cnt) = NMFCluster(H{i}{cnum});
               cnt = cnt + 1;
           end
       end
   end   
   
   % Prepare empty consensous matrices 
    for cnum = kvect
        if (any(strcmp(cc_type,'hard')))
            cc_clust.clust_hard.count_co_clust{cnum} = nan(Nind);
        end
        if (any(strcmp(cc_type,'soft')))            
            cc_clust.clust_soft.count_co_clust{cnum} = nan(Nind);
        end
        if (any(strcmp(cc_type,'euclid')))
            cc_clust.clust_soft_euclid.count_co_clust{cnum} = nan(Nind);
        end
        if (any(strcmp(cc_type,'corr')))
            cc_clust.clust_soft_corr.count_co_clust{cnum} = nan(Nind);
        end
        if (any(strcmp(cc_type,'kldiv')))
            cc_clust.clust_soft_kldiv.count_co_clust{cnum} = nan(Nind);
        end
    end    
    cc_clust.cc_count_all = zeros(Nind);    

    % Summarize results of clustering
    for i = 1:nsample
        cnt = 0;
        for cnum = kvect
            cnt = cnt + 1;
            
            hard_label = hard_label_vect{i}(:,cnt);
            Nsample = Nsample_vect{i};
            cH = H{i}{cnum};
            
            if (~isempty(outlierVect))
                hard_label = hard_label(outlierVect{i}{cnum});
                Nsample = Nsample(outlierVect{i}{cnum});
                cH = cH(:,outlierVect{i}{cnum});
            end
            
            if (any(strcmp(cc_type,'hard')))
                % nanadd(,cc_all.clust_hard.count_co_clust{z},cc_cluster.clust_hard.count_co_clust{z}),3);
                cc_clust.clust_hard.count_co_clust{cnum} = ...
                    nanadd(cc_clust.clust_hard.count_co_clust{cnum},...
                    label_to_ccmat_quick(hard_label,Nsample,Nind));
            end
            if (any(strcmp(cc_type,'soft')))                
                cc_clust.clust_soft.count_co_clust{cnum} = ...
                    nanadd(cc_clust.clust_soft.count_co_clust{cnum},...
                    label_to_ccmat_soft(cH,Nsample,Nind));
            end
            if (any(strcmp(cc_type,'euclid')))
                cc_clust.clust_soft_euclid.count_co_clust{cnum} = ...
                    nanadd(cc_clust.clust_soft_euclid.count_co_clust{cnum},...
                    label_to_ccmat_soft_euclid(cH,Nsample,Nind));
            end
            if (any(strcmp(cc_type,'corr')))
                cc_clust.clust_soft_corr.count_co_clust{cnum} = ...
                    nanadd(cc_clust.clust_soft_corr.count_co_clust{cnum},...
                    label_to_ccmat_soft_corr(cH,Nsample,Nind));
            end
            if (any(strcmp(cc_type,'kldiv')))
                cc_clust.clust_soft_kldiv.count_co_clust{cnum} = ...
                    nanadd(cc_clust.clust_soft_kldiv.count_co_clust{cnum},...
                    label_to_ccmat_soft_kldiv(cH,Nsample,Nind));
            end
            
        end
        cc_clust.cc_count_all = cc_clust.cc_count_all + perm_to_count(Nsample,Nind);
    end
    
    for cnum = kvect
        if (any(strcmp(cc_type,'hard')))
            cc_clust.clust_hard.cc_matrix{cnum} = (cc_clust.clust_hard.count_co_clust{cnum})./(max(cc_clust.cc_count_all,1));
        end
        if (any(strcmp(cc_type,'soft')))
            cc_clust.clust_soft.cc_matrix{cnum} = (cc_clust.clust_soft.count_co_clust{cnum})./(max(cc_clust.cc_count_all,1));
        end
        if (any(strcmp(cc_type,'euclid')))
            cc_clust.clust_soft_euclid.cc_matrix{cnum} = (cc_clust.clust_soft_euclid.count_co_clust{cnum})./(max(cc_clust.cc_count_all,1));
        end
        if (any(strcmp(cc_type,'corr')))
            cc_clust.clust_soft_corr.cc_matrix{cnum} = (cc_clust.clust_soft_corr.count_co_clust{cnum})./(max(cc_clust.cc_count_all,1));
        end
        if (any(strcmp(cc_type,'kldiv')))
            cc_clust.clust_soft_kldiv.cc_matrix{cnum} = (cc_clust.clust_soft_kldiv.count_co_clust{cnum})./(max(cc_clust.cc_count_all,1));
        end
    end    
end