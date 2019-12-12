function [matrix] = convertCellArrayToMat(cellArray)
    sizeOfData = size(cellArray);
    matrix = zeros(sizeOfData(1),sizeOfData(2));
    for k = 1:sizeOfData(1)
        for j = 1:sizeOfData(2)
            matrix(k, j) = cellArray{k,j};
        end
    end
end