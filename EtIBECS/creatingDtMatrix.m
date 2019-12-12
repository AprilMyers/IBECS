function [dtMatrix] = creatingDtMatrix(Matrix)
    sizeOfData = size(Matrix);
    dtMatrix = zeros(sizeOfData(1), sizeOfData(2) - 1);
    
    for trialNumber = 1:sizeOfData(1)
        for frameNumber = 1:(sizeOfData(2) - 1)
            dtMatrix(trialNumber, frameNumber) = Matrix(trialNumber, frameNumber) - Matrix(trialNumber, frameNumber + 1);
        end
    end
end

