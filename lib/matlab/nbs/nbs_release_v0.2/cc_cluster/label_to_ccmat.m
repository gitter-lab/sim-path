% Count co-clustering pairs
function co_cluster = label_to_ccmat(tlabel_exp,perm_list,Nind)
    co_cluster = nan(Nind);
    
    try     
        [sperm,sIdx] = sort(perm_list);
        
        for i = 1:length(sperm)
            zi = sperm(i);
            ziIdx = sIdx(i);
            co_cluster(zi,zi)  = 1;
            
            for j=(i+1):length(sperm)
                co_cluster(zi,sperm(j)) = tlabel_exp(ziIdx) == tlabel_exp(sIdx(j));
            end
        end
    catch e
        fprintf(1,'Warn: bad_cluster assignment\n');
    end
    
    co_cluster = nanmax(co_cluster,co_cluster');
end
