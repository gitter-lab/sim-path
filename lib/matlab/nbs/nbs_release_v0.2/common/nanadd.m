function sumM = nanadd(A,B)
% 
% NaN ignoring add. 
% The sum of x + nan = x and nan + nan = nan  
%
% 

    zsum = cat(3,A,B);
    znan = all(isnan(zsum),3);
    
    zsum(isnan(zsum)) = 0;
    sumM = sum(zsum,3);
    sumM(znan) = nan;
    
end

