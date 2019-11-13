function [minX,maxX,minY,maxY] = selectCroppedRegion(cFNames, iTOI, dFldr)
%SELECTCROPPEDREGION Summary of this function goes here
%   Detailed explanation goes here

    pupilFile = strtrim([dFldr cFNames{iTOI}]);
    fullImage = imread(pupilFile, 1);
    figure
    
    acceptCrop = 0;
    while acceptCrop==0
        beep
        [J,rect2] = imcrop(fullImage*4), hold on

        rect2 = round(rect2);
        minX = rect2(1);
        maxX = rect2(1)+rect2(3);
        minY = rect2(2);
        maxY = rect2(2)+rect2(4);

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