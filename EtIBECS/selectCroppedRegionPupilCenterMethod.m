function [minX,maxX,minY,maxY] = selectCroppedRegionPupilCenterMethod(cFNames, iTOI, dFldr)
%SELECTCROPPEDREGION Summary of this function goes here
%   Detailed explanation goes here

    pupilFile = strtrim([dFldr cFNames{iTOI}]);
    fullImage = imread(pupilFile, 1);
    
    
    acceptCrop = 0;
    while acceptCrop==0
        beep
        imshow(fullImage);
        [X, Y] = ginput(1), hold on

        X = round(X);
        Y = round(Y);
        
        minX = X - 50;
        maxX = X + 50;
        minY = Y - 50;
        maxY = Y + 50;

        % scatter([minX minX maxX maxX], [maxY minY minY maxY],'r');
        xs = [minX maxX maxX minX minX];
        ys = [minY minY maxY maxY minY];
        figure
        plot(xs, ys,'r','linestyle','-');
        cropImage = fullImage([minY:maxY],[minX:maxX]);
        title('full')
        imagesc(cropImage)
        % caxis([0 5])
        % colormap hot
        title('cropped')
        beep
        acceptCrop = input('Accept Crop? Enter 1 to accept or 0 to reject: ')
    end
end