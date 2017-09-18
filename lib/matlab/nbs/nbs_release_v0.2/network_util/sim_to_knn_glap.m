function [knnGlap,knnMat] = sim_to_knn_glap(simMatrix,kn,isDist)

    ordSimOrDist = 'descend';
    if ( exist('isDist','var') ~= 0  && isDist == 1)
        ordSimOrDist = 'ascend';
    end
    
    knnMat = false(size(simMatrix));
    
    [C,sIdx] = sort(simMatrix,2,ordSimOrDist);
    colMat = repmat((1:length(simMatrix))',1,kn);
    sIdxSelect = sIdx(:,1:kn);;
    
    
    linInd = sub2ind(size(simMatrix),sIdxSelect(:),colMat(:));
    
    knnMat(linInd) = 1;
    
    knnMat = or(knnMat,knnMat');
    
    knnGlap = diag(sum(knnMat)) - knnMat;
   

end