function pupil(inputArg1,inputArg2)
%PUPIL Summary of this function goes here
%   Detailed explanation goes here

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
    [minX, maxX, minY, maxY] = selectCroppedRegion(cFNames, iTOI);
        
    % Call Pupil Analysis function on each trial
    tic
    parfor u = 1:nTrls
        [pupilSizeTrialU, pupilVids] = pupilAnalysis(dFldr, u, iTOI, cFNames, minX, maxX, minY, maxY);
        pupilSizes{u} = pupilSizeTrialU;
        % Optional to draw figures
        visualizePupilAnalysis(minX, maxX, minY, maxY);
    end
    toc
    
    % Converts the cell array to a Matrix of all values and use Nans to pad
    [pupilSizesMat] = convertCellArraytoMat(pupilSizes);
    
    % Create the dt matrix
    [pupilSizedt] = creating_dt(pupilSizes);
    
    % Filter the dt array
    [pupilareaProc, pupilsize_filtered] = filterMatrix(pupilSizes, pupilSizedt);
    
    % Converts the cell array to a matrix of all values and use Nans to pad
    [pupilsize_filtered_mat] = convertCellArraytoMat(pupilsize_filtered);
    
    % Interpolation
    [pupilsize_filtered] = interpolation(pupilsize_filtered);
    
     % Converts the cell array to a matrix of all values and use Nans to pad
    [pupilsize_interp_mat] = convertCellArraytoMat(pupilsize_filtered);
    
end

