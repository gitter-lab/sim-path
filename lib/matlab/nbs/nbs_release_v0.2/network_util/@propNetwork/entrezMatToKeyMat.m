function [fullKeyMat,unMappedGenoMat] = entrezMatToKeyMat(netobj,igeno,genoEntrezKey)
    
    
    fullKeyMat = zeros(size(igeno,1),length(netobj.key));
    unMappedGenoMat = zeros(size(igeno));
    for i = 1:size(igeno,1)
        nnzPos = logical(igeno(i,:));
        nnzPosIdx = find(nnzPos);
        [tt,tunmapped] = netobj.mapEntrezListToKeyVect(genoEntrezKey(nnzPos),igeno(i,nnzPos));
        
        fullKeyMat(i,:) = tt;    
        unMappedGenoMat(i,nnzPosIdx(tunmapped)) = 1;
    end


end