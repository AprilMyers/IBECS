function [minX,maxX,minY,maxY] = selectCroppedRegionPupilCenterMethod(tifFiles, dFldr)
%SELECTCROPPEDREGION Summary of this function goes here
%   Detailed explanation goes here

    pupilFile = strtrim([dFldr tifFiles{1}]);
    disp(tifFiles{1});
    fullImage = imread(pupilFile, 1);
    
    
    acceptCrop = 0;
    while acceptCrop==0
        beep
        imshow(fullImage);
        [X, Y] = ginputWhite(1);
        hold on;

        X = round(X);
        Y = round(Y);
        
        minX = X - 40;
        maxX = X + 40;
        minY = Y - 40;
        maxY = Y + 40;

        % scatter([minX minX maxX maxX], [maxY minY minY maxY],'r');
        xs = [minX maxX maxX minX minX];
        ys = [minY minY maxY maxY minY];
        figure
        plot(xs, ys,'r','linestyle','-');
        cropImage = fullImage([minY:maxY],[minX:maxX]);
        title('full');
        imagesc(cropImage);
        % caxis([0 5])
        % colormap hot
        title('cropped');
        beep
        acceptCrop = input('Accept Crop? Enter 1 to accept or 0 to reject: ');
    end
    close all;
end