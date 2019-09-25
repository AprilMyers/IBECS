function [ contStr ] = cellContainsStr( A, str )
% finds cells within array that contains the string
contStr = ~cellfun(@isempty,strfind(A,str));
end

