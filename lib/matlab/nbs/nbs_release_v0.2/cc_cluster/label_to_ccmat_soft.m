% Label to cc mat soft
function [co_cluster] = label_to_ccmat_soft(H,perm_list,Nind)
    co_cluster = nan(Nind);
    
    Hnorm = bsxfun(@times,H,max(sum(H),eps).^-1);
    Hnorm(:,all(isnan(Hnorm)))=1/size(H,1);
    
    [sperm,sIdx] = sort(perm_list);
    cc_weight = Hnorm'*Hnorm;

    co_cluster(sperm,sperm) = cc_weight(sIdx,sIdx);
end
