function [keyGeneVect,unmappedGeneVect] = mapEntrezListToKeyVect(netobj,inEntrezList,inValueList)
%
% Map entrez list to Key vect
% 
    if (exist('inValueList','var') == 0)
        inValueList = ones(size(inEntrezList));
    end
    % Map to internal key
    keyList = nanvalues(netobj.entrezIDtoKeyMap,inEntrezList');
        
    keyGeneVect = zeros(1,length(netobj.key));
    unmappedGeneVect = false(1,length(inEntrezList));
    
    if (iscell(keyList))
        for i = 1:length(keyList)
            if (isempty(keyList{i}))                
                unmappedGeneVect(i) = 1;
                continue;
            end
            
            keyListP = keyList(i);
            
            kPos = nanvalues(netobj.keyPosMap,keyListP);
            if (isnan(kPos))
                unmappedGeneVect(i) = 1;
                fprintf(1,'Why is this?\n');
                continue;
            end
            
            if (length(kPos) == 1)
                keyGeneVect(kPos(1)) =  keyGeneVect(kPos(1)) + inValueList(i);
            else
                warn('Ambigious ID');
                keyGeneVect(kPos{1}{1}) =  keyGeneVect(kPos{1}{1}) + inValueList(i);
%                 for zpos = kPos
%                     keyGeneVect(zpos{:}) = keyGeneVect(zpos{:}) + 1/length(kPos);
%                 end
            end            
        end
    else   
        kPos = nanvalues(netobj.keyPosMap,keyList);

        kPosNan = isnan(kPos);
        
        if (all(kPosNan))                
            unmappedGeneVect(:) = 1;
            fprintf(1,'Why is this?\n');            
        else
            keyGeneVect(kPos(~kPosNan)) = 1;
            unmappedGeneVect(kPosNan) = 1;
        end
        
%         for i = 1:length(keyList)
%             if (isempty(keyList(i)) || isnan(keyList(i)) )
%                 unmappedGeneVect(i) = 1;                
%                 continue;
%             end
%             
%             keyListP = keyList(i);
%             
%             kPos = nanvalues(netobj.keyPosMap,keyListP);
%             
%             if (isnan(kPos))                
%                 unmappedGeneVect(i) = 1;
%                 fprintf(1,'Why is this?\n');
%                 continue;
%             end
%             
%             if (length(kPos) == 1)
%                 keyGeneVect(kPos) = keyGeneVect(kPos) + inValueList(i);
%             else
%                 warn('Ambigious ID');
%                 keyGeneVect(kPos(1)) = keyGeneVect(kPos(1)) + inValueList(i);
% %                 for zi = 1:length(kPos)
% %                     zpos = kPos(zi);
% %                     keyGeneVect(zpos) = keyGeneVect(zpos) + 1/length(kPos);
% %                 end
%             end            
%         end
    end
    
    

end 