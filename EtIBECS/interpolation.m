function [pupilsize_filtered] = interpolation(pupilsize_filtered)
    %nterpolate through pupilsize_filtered

    % Outer for loop to iterate between different videos
    for u = 1:length(pupilsize_filtered)


        %Inner for loop to iterate between different frames seraching for a Nan
        for k = 1:length(pupilsize_filtered{1, u})
            %counters to increment up or down frames to find the nearrest
            %values in the while loops
            j = 0;
            m = 0;
            %Variables to save the nearest values in
            closest_prev = 0;
            closest_next = 0;

            if isnan(pupilsize_filtered{1, u}(k))

                %find previous frame area 
                while k - j > 0
                    if (k - j) == 1
                        closest_prev = NaN;
                        break;
                    end
                    if ~isnan(pupilsize_filtered{1, u}(k - j))
                        closest_prev = pupilsize_filtered{1, u}(k-j);
                    end
                    j = j + 1;
                end

                %find next frame area
                while k + m < (length(pupilsize_filtered{1, u}) + 1)
                    if (k + m) == length(pupilsize_filtered{1, u})
                        closest_next = NaN;
                        break;
                    end
                    if ~isnan(pupilsize_filtered{1, u}(k + m))
                        closest_next = pupilsize_filtered{1, u}(k+m);
                        break;
                    end
                    m = m + 1;
                end

                %average closest_next and closest_prev to get value for
                %previous nan valued frame
                if isnan(closest_prev)
                    average = closest_next;
                elseif isnan(closesy_next)
                    average = closest_prev;
                else
                    average = closest_prev + closest_next;
                end
                pupilsize_filtered{1, u}(k) = average;
            end
        end

    end
end