function [nTrls, dFldr, cFNames, aTrl, sTrl, nCam1, nCam2, nLick, nMWhl,... 
    nTrig, nValv, nShck, nMove, nProt, nStim, fSesh, ... 
    fCam1, fCam2, fLick, fMWhl, fTrig, fValv, fShck, ...
    fMove, fProt, fStim, fTxt, fTiff, fSTrl, seshTitle] = setUpPupilAnalysis()
%SETUPPUPILANALYSIS Summary of this function goes here
%   Detailed explanation goes here

    % Get the directory with the videos
    dFldr = uigetdir();
    
    %Session name
    seshN = 's';  
    
    % Sorting the Files in the Directory
    % Trial types
    aTrl = 'adpt';
    sTrl = 'stim';
    
    % Measurement types
    nCam1 = 'Camera 1';
    nCam2 = 'Camera 2';
    nLick = 'Lick Port';
    nMWhl = 'Mouse Wheel';
    nTrig = 'Trigger';
    nValv = 'Valve';
    nShck = 'Shock';
    nMove = 'Move';
    nProt = 'Protocol';
    nStim = 'angle';
    
    % find files
    dFldr = [dFldr '\'];
    cFNames = num2cell(lswindows(dFldr),2); %cell array of file names
    fSesh = cellContainsStr(cFNames,seshN); % Session files
    fCam1 = cellContainsStr(cFNames,nCam1); %Camera 1 files
    fCam2 = cellContainsStr(cFNames,nCam2);
    fLick = cellContainsStr(cFNames,nLick);
    fMWhl = cellContainsStr(cFNames,nMWhl);
    fTrig = cellContainsStr(cFNames,nTrig);
    fValv = cellContainsStr(cFNames,nValv);
    fShck = cellContainsStr(cFNames,nShck);
    fMove = cellContainsStr(cFNames,nMove);
    fProt = cellContainsStr(cFNames,nProt);
    fStim = cellContainsStr(cFNames,nStim);
    fTxt = cellContainsStr(cFNames,'txt'); % file type
    fTiff = cellContainsStr(cFNames,'tif');
    
    fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',1)]);
    
    %seshTitle = cFNames{10}(1:find('_'==cFNames{10},1,'first')-1);
    seshTitle = "This is the seshTitle";
    
    iTOI = fProt;
    fNm = strtrim([dFldr cFNames{iTOI}]);
    protocol = readtable(fNm, 'ReadVariableNames', true, 'Delimiter', '\t');
    
    iTOI = find(fStim);
    fNm = strtrim([dFldr cFNames{iTOI}]);
    stimInfo = load(fNm);
    
    % number of trials
%     nTrls = str2num(protocol.TrialName{1});
     % Calculating nTrls using for loop to count the number of tif files
     
%      nTrls = sum(fTiff);
%      disp(["SETUP's value for nTrls " nTrls])
    nTrls = sum(fTiff);
    nTrlsCheck = length(protocol.TrialName)-1;
    if nTrls~=nTrlsCheck
        warning('Mismatched number of trials in Protocol')
    end
end