function  [new_adj_mat,gvect] = extractThrehold(netobj,quantile_threshold)
    
    if (netobj.qthreshold ~= 0)
        error('Network already thresheld');
    end

    sval_list = sort(nonzeros(netobj.adj_mat));
    q_threshold = sval_list(floor(length(sval_list)*quantile_threshold));
    
    new_adj_mat = netobj.adj_mat;    
    new_adj_mat(netobj.adj_mat < q_threshold) = 0;
       
    tt = all(new_adj_mat == 0);
    new_adj_mat = new_adj_mat(~tt,~tt);
    gvect = netobj.key(~tt);
    
%     [expmap_keys,expmap_values] = expanFullMap(netobj.KeytoEntrezIDMap); 
%     if (iscell(expmap_values))
%         expmap_values = cell2mat(expmap_values);
%     end
%         
%     [listSelect,idxlist] = selectFromList(expmap_keys,gvect);
%     
%     entrez2key = createFullMap(expmap_values(idxlist),listSelect);
%     key2entrez = createFullMap(listSelect,expmap_values(idxlist));
%     
%     obj = multiIDNetwork(new_adj_mat,gvect,entrez2key,key2entrez);
%     obj.qthreshold = quantile_threshold;
end 