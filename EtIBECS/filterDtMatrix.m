function [dtMatrixFiltered] = filterDtMatrix(dtMatrix, threshold)
    sizeOfData = size(dtMatrix);
    dtMatrixFiltered = zeros(sizeOfData(1), sizeOfData(2));
    for trialNumber = 1:sizeOfData(1)
        for frameNumber = 1:sizeOfData(2)
            if dtMatrix(trialNumber, frameNumber) < threshold && dtMatrix(trialNumber, frameNumber) > -threshold
                dtMatrixFiltered(trialNumber, frameNumber) = dtMatrix(trialNumber, frameNumber);
            else
                dtMatrixFiltered(trialNumber, frameNumber) = nan;
            end
        end
    end
end
