% Count co-clustering pairs
function co_cluster = label_to_ccmat_quick(tlabel_exp,perm_list,Nind)
    co_cluster = nan(Nind);
    tlabel_full = nan(Nind,1);
    
    [sperm,sIdx] = sort(perm_list);    
    tlabel_full(sperm) = tlabel_exp(sIdx);
    
    label_list = unique(tlabel_full);
    label_list = label_list(~isnan(label_list));
    label_list = label_list(:)';
    
    znan = ~isnan(tlabel_full);
    co_cluster(znan,znan) = 0;
    for ck = label_list
        co_cluster(tlabel_full == ck,tlabel_full == ck) = 1;
    end
    
end
