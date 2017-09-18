classdef propNetwork 
    properties
        adj_mat
        adj_mat_norm = [];
        key
        keyPosMap
        qthreshold = 0;
        
        entrezIDtoKeyMap
        symbolIDtoKeyMap
        KeytoEntrezIDMap
       
        adj_mat_norm_val = [];
        propVal = [];
    end
    
    methods(Static)
        normD = normAdjMat(adj_mat);
        [listSelect,idxlist] = selectFromList(listAll,listKey)
    end
    methods                                
        function obj = propNetwork(adj_mat,key,entrezIDtoKeyMap,KeytoEntrezIDMap,symbolIDtoKeyMap)
            obj.adj_mat = adj_mat;
            obj.key = key;
            obj.entrezIDtoKeyMap = entrezIDtoKeyMap;
            obj.symbolIDtoKeyMap = symbolIDtoKeyMap;
            obj.KeytoEntrezIDMap = KeytoEntrezIDMap;
            obj.keyPosMap = containers.Map(key,1:length(key));
            
            % obj.adj_mat_norm = propNetwork.normAdjMat(adj_mat);
        end         
        
        [adj_mat,gvect] = extractThrehold(netobj,alpha)
        
        [keyGeneVect,unmappedGeneVect] = mapEntrezListToKeyVect(netobj,inEntrezList,inValuesList)        
        
        inEntrezList = mapKeyVectToEntrezList(netobj,keyGeneVect)
        
        [expgeno,cover_fraction,cid] = networkPropgateEntrez(netobj,genoEntrez,genoEntrezKey,alpha,isInclude,doNorm,extKMat)
        
        [fullKeyMat,unMappedGenoMat] = entrezMatToKeyMat(netobj,igeno,genoEntrezKey)
        
        [expgeno,cover_fraction,cid] = networkPropgateEntrezQuick(netobj,genoEntrez,genoEntrezKey,alpha,isInclude,doNorm)
        
        obj = generateQuickPropmat(netobj,netpropFactor)
    end
    

end