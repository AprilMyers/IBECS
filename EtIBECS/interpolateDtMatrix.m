function [interpolatedMatrix] = interpolateDtMatrix(inputMatrix)
    %nterpolate through pupilsize_filtered
    sizeOfData = size(inputMatrix);
    interpolatedMatrix = zeros(sizeOfData(1), sizeOfData(2));
    
    for trialNumber = 1:sizeOfData(1)
        for frameNumber = 1:sizeOfData(2)
            j = 0;
            m = 0;
            
            closest_prev = 0;
            closest_next = 0;
            
            
            if isnan(inputMatrix(trialNumber, frameNumber))
                
                if frameNumber == 1
                    closest_prev = NaN;
                else %~isnan(inputMatrix(trialNumber, frameNumber - 1))
                    closest_prev = inputMatrix(trialNumber, frameNumber - 1);
                end
                    
                k = 0;
                closest_next_index = frameNumber - 1;
                while frameNumber + m < sizeOfData(2) + 1 && k == 0
                    if (frameNumber + m) == sizeOfData(2)
                        closest_next = NaN;
                        k = 1;
                    elseif ~isnan(inputMatrix(trialNumber, frameNumber + m))
                        closest_next = inputMatrix(trialNumber, frameNumber + m);
                        
                        closest_next_index = frameNumber + m;
                        k = 1;
                    end
                    m = m + 1;
                end
                
                
                if isnan(closest_prev)
                    average = closest_next;
                elseif isnan(closest_next)
                    average = closest_prev;
                elseif isnan(closest_prev) && isnan(closest_next)
                    average = 0;
                else
                    average = (closest_prev + closest_next)/2;
                end
                for d = frameNumber:closest_next_index - 1
                    interpolatedMatrix(trialNumber, d) = average;
                    inputMatrix(trialNumber, d) = average;
                end
            else
                interpolatedMatrix(trialNumber, frameNumber) = inputMatrix(trialNumber, frameNumber);
            end
        end
    end
end