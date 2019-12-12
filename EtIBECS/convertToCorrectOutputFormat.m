function [data] = convertToCorrectOutputFormat(dataForAllTrialsArray, minFrames, nTrls)
    Areas = cell(nTrls, minFrames);
    Centers = cell(nTrls, minFrames);
    MajorAxisLength = cell(nTrls, minFrames);
    MinorAxisLength = cell(nTrls, minFrames);
    BoundingBox = cell(nTrls, minFrames);
    
    matrixSize = size(Areas);
    for k = 1:matrixSize(1) %number of trials
        for j = 1:matrixSize(2) %number of frames
            trialData = dataForAllTrialsArray{k,1};
            frameData = trialData{j,1};
            frameArea = frameData.Area;
            frameCenter = frameData.Centroid;
            frameMajorAxisLength = frameData.MajorAxisLength;
            frameMinorAxisLength = frameData.MinorAxisLength;
            frameBoundingBox = frameData.BoundingBox;
            
            Areas(k,j) = {frameArea};
            Centers(k,j) = {frameCenter};
            MajorAxisLength(k,j) = {frameMajorAxisLength};
            MinorAxisLength(k,j) = {frameMinorAxisLength};
            BoundingBox(k,j) = {frameBoundingBox};
        end
    end
    
    data = struct;
    data.AreasArray = Areas;
    data.CentersArray = Centers;
    data.MajorAxisLengthArray = MajorAxisLength;
    data.MinorAxisLengthArray = MinorAxisLength;
    data.BoundingBoxArray = BoundingBox;
end