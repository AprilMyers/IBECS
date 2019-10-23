function [contStr] = cellContainsStr(A, str)
%finds the cells in a cell array that contain the string str
    contStr = ~cellfun(@isempty,strfind(A,str));
end