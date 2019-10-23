function [pupilSizedt] = creating_dt(pupilSizes)
    % Generate dt list
    pupilSizedt = cell(1,91);
    for u = 1:91 
        pupilsize = pupilSizes{u};
        t = 1:length(pupilsize);
        v = zeros(length(t)-1,1);
        for i = 1:length(t)-1
          v(i) = (pupilsize(i+1)-pupilsize(i))/(t(i+1)-t(i));
          pupilSizedt{u} = v ;
        end
    end

    %This is cell array of dt values
    %pupilSizedt
end