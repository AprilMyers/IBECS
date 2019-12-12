function [cropImage, data] = pupil()
%PUPIL: Outputs a structure containing 4 fields
% 1) Area
% 2) MajorAxisLength 
% 3) xPosition
% 4) yPosition
    % Decide whether you want to create plots
    
    
    createPlot = input("Do you want to create plots? 1 to create, 0 to not create: ");
    % Decide whether you want to display plots
    displayPlot = input("Would you like to display plots? 1 to display, 0 to not display: ");
    
    if displayPlot == 1
        howManyPlot = input("How many plots would you like to display? (all is not recommended) ");
        if howManyPlot == 0
            howManyPlot = 1;
        end
    else
        howManyPlot = 1;
        displayPlot = "off";
    end
    
    % Call Setup and File Read In Function
    disp("Started Setup");
    [nTrls, dFldr, cFNames, aTrl, sTrl, nCam1, nCam2, nLick, nMWhl,... 
        nTrig, nValv, nShck, nMove, nProt, nStim, fSesh, ... 
        fCam1, fCam2, fLick, fMWhl, fTrig, fValv, fShck, ...
        fMove, fProt, fStim, fTxt, fTiff, fSTrl,seshTitle] = setUpPupilAnalysis();
    
    
    
    disp("Finished Setup");
    disp("Hello")
    disp(dFldr)
    if ~exist(fullfile(dFldr, 'Figures'),'dir')
        mkdir(fullfile(dFldr, 'Figures'))
    end
    
    %Ryans Figure Settings to Display Figures
    RyansFigureSettings;
    
    
    % Filter files based on name attributes, tifFiles holds the tif files used to generate the pupil analysis
    iTOI = find(fSesh & fCam1 & fSTrl & fTiff);
    
    tifFiles = cell(sum(fTiff), 1); 
    count = 1;
    for k = 1:size(fTiff)
        if fTiff(k) == 1
            tifFiles{count} = cFNames{k};
            count = count + 1;
        end
    end
%     disp(cFNames);
    disp("hello");
    disp(iTOI);
%     disp(tifFiles);
        
    % Call Cropping function
    disp("Started Cropping");
    [minX, maxX, minY, maxY] = selectCroppedRegionPupilCenterMethod(tifFiles, dFldr);
    disp("Finished Cropping");
    % Call Pupil Analysis function on each trial
    tic
    
    
    %for debugging%
    nTrls = 1;
    %end debugging%
    
    dataForAllTrialsArray = cell(nTrls, 1);
    minFrames = 10000000;
    for trialNumber = 1:nTrls
        disp("Started pupilAnalysis for trial:");
        disp(trialNumber);
        
%         [dataTrialArray] = pupilAnalysis(dFldr, tifFiles, trialNumber, minX, maxX, minY, maxY);
        [cropImage, dataTrialArray, nFrames] = pupilAnalysis(howManyPlot, displayPlot, createPlot, dFldr, tifFiles, trialNumber, minX, maxX, minY, maxY);
        dataForAllTrialsArray(trialNumber, 1) = {dataTrialArray};
        
        if nFrames < minFrames
            minFrames = nFrames;
        end
        
        disp("Done with trial:");
        disp(trialNumber);
        disp("Corresponding file:");
        disp(strtrim(tifFiles{trialNumber}));
    end
    toc
    
% Convert dataTrialArray into struct with the following fields:
    % --> Areas (Trial x Frame Matrix of Areas)
    % --> Centers (Trial X Frame Matrix of Centers)
    % --> MajorAxisLength (Trial x Frame Matrix of MajorAxisLength)
    % --> MinorAxisLength (Trial X Frame Matrix of MinorAxisLength)
    % --> BoundingBox (Trial x Frame Matrix of BoundingBox)
    
    disp("Started converting to correct output format");
    data = convertToCorrectOutputFormat(dataForAllTrialsArray, minFrames, nTrls);
   
    disp("Started converting from cell array of areas to matrix of areas");
    data.AreasMatrix = convertCellArrayToMat(data.AreasArray);
    
    disp("Started calculating dtAreas matrix");
    data.dtAreasMatrix = creatingDtMatrix(data.AreasMatrix);
    
    disp("Started filtering the dtAreas matrix");
    threshold = std(data.dtAreasMatrix); %Needs to be adjusted, can be adjusted later
    data.dtAreasMatrixFiltered = filterDtMatrix(data.dtAreasMatrix, threshold);
    
    disp("Started calculating dtAreasInterpolated matrix");
    data.dtAreasMatrixFilteredInterpolated = interpolateDtMatrix(data.dtAreasMatrixFiltered);
    
    save('all_data.mat','-struct','data') %saves to pupilAnalysisNew folder

end