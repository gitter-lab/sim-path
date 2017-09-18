% Count all pairs
function count_mat = perm_to_count(perm_N,Nind)
  
    count_mat = zeros(Nind);
        
    count_mat(perm_N,perm_N) = 1;
         
end
