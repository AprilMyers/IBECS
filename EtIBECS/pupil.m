% function pupil(inputArg1,inputArg2)
% %PUPIL Summary of this function goes here
% %   Detailed explanation goes here
% 
%     % Call Setup and File Read In Function
%     [nTrls, dFldr, cFNames, aTrl, sTrl, nCam1, nCam2, nLick, nMWhl,... 
%         nTrig, nValv, nShck, nMove, nProt, nStim, fSesh, ... 
%         fCam1, fCam2, fLick, fMWhl, fTrig, fValv, fShck, ...
%         fMove, fProt, fStim, fTxt, fTiff, fSTrl,seshTitle] = setUpPupilAnalysis();
%     
%     %Ryans Figure Settings to Display Figures
%     RyansFigureSettings;
%     
%     
%     % Filter files based on name attributes
%     iTOI = find(fSesh & fCam1 & fSTrl & fTiff);
%         
%     % Call Cropping function
%     [minX, maxX, minY, maxY] = selectCroppedRegionPupilCenterMethod(cFNames, iTOI, dFldr);
%         
%     % Call Pupil Analysis function on each trial
%     tic
%     pupilSizes = {};
%     for u = 1:1 %1:nTrls
%         [pupilSizeTrialU, pupilVids] = pupilAnalysis(dFldr, u, iTOI, cFNames, minX, maxX, minY, maxY);
%         pupilSizes{u} = pupilSizeTrialU;
%         disp("Done with trial");
%         disp(u);
%     end
%     toc
%     
%     % Converts the cell array to a Matrix of all values and use Nans to pad
%     disp("The dimension of pupilSizes is: ")
%     pupilSizes
%     
% %     [pupilSizesMat] = convertCellArraytoMat(pupilSizes);
% %     
% %     % Create the dt matrix
% %     [pupilSizedt] = creating_dt(pupilSizes);
% %     
% %     % Filter the dt array
% %     [pupilareaProc, pupilsize_filtered] = filterMatrix(pupilSizes, pupilSizedt);
% %     
% %     % Converts the cell array to a matrix of all values and use Nans to pad
% %     [pupilsize_filtered_mat] = convertCellArraytoMat(pupilsize_filtered);
% %     
% %     % Interpolation
% %     [pupilsize_filtered] = interpolation(pupilsize_filtered);
% %     
% %      % Converts the cell array to a matrix of all values and use Nans to pad
% %     [pupilsize_interp_mat] = convertCellArraytoMat(pupilsize_filtered);
%     
% end



function [cFNames iTOI, pupilS] = pupil(inputArg1,inputArg2)
%PUPIL: Outputs a structure containing 4 fields
% 1) Area
% 2) MajorAxisLength 
% 3) xPosition
% 4) yPosition


    % Call Setup and File Read In Function
    [nTrls, dFldr, cFNames, aTrl, sTrl, nCam1, nCam2, nLick, nMWhl,... 
        nTrig, nValv, nShck, nMove, nProt, nStim, fSesh, ... 
        fCam1, fCam2, fLick, fMWhl, fTrig, fValv, fShck, ...
        fMove, fProt, fStim, fTxt, fTiff, fSTrl,seshTitle] = setUpPupilAnalysis();
    
    %Ryans Figure Settings to Display Figures
    RyansFigureSettings;
    
    
    % Filter files based on name attributes
    iTOI = find(fSesh & fCam1 & fSTrl & fTiff);
        
    % Call Cropping function
    [minX, maxX, minY, maxY] = selectCroppedRegionPupilCenterMethod(cFNames, iTOI, dFldr);
        
    % Call Pupil Analysis function on each trial
    tic
    disp("This is the number of trials: ")
    disp(nTrls)
    mkdir(dFldr, "Figures");
    for u = 1:1 %nTrls
        [pupilSizeTrialU, pupilVids] = pupilAnalysis(dFldr, u, iTOI, cFNames, minX, maxX, minY, maxY);
        
%         disp("Trial rpData_pam ", u);
%         rpData_pam
%         
%         
%         disp("Trial rpData_wholething ", u);
%         rpData_wholething
%         
        pupilSizes{u} = pupilSizeTrialU;
        disp("Done with trial");
        disp(u);
    end
    toc
    
    pupilS = pupilSizes;
%     % Converts the cell array to a Matrix of all values and use Nans to pad
%     [pupilSizesMat] = convertCellArraytoMat(pupilSizes);
%     
%     % Create the dt matrix
%     [pupilSizedt] = creating_dt(pupilSizes);
%     
%     % Filter the dt array
%     [pupilareaProc, pupilsize_filtered] = filterMatrix(pupilSizes, pupilSizedt);
%     
%     % Converts the cell array to a matrix of all values and use Nans to pad
%     [pupilsize_filtered_mat] = convertCellArraytoMat(pupilsize_filtered);
%     
%     % Interpolation
%     [pupilsize_filtered] = interpolation(pupilsize_filtered);
%     
%      % Converts the cell array to a matrix of all values and use Nans to pad
%     [pupilsize_interp_mat] = convertCellArraytoMat(pupilsize_filtered);
    
end