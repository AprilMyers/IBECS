function [pupilsize_filtered, pupilareaProc] = filterMatrix(pupilSizes, pupilSizedt)
    pupilsize_filtered = cell(1,91);

    for index=1:91
        disp(index)
        pupilarea_dt = single(pupilSizedt{index}) ;
        pupilarea_dt = [NaN;pupilSizedt{index}(1:end)];

        area_ses = pupilSizes{index};
        %ses1 = table(area_ses1_dt, area_ses1);
        % extract first row: ses1(1:501, 1)
        isgoodframe = (-1.5 < pupilarea_dt & pupilarea_dt < 1.5) ;
        area_ses(~isgoodframe) = nan;

        pupilsize_filtered{index} = area_ses;
    end

    pupilareaProc = pupilsize_filtered;  
    %%  

    % Plotting stuff save for later
%     num_sesplot = 1
%     plot(pupilSizes{num_sesplot})
%     hold on
%     plot(pupilsize_filtered{num_sesplot})
%     hold off
end
