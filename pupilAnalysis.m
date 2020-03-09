function [o2, p2, r2, c2, CH2, cropImage, dataTrialArray, nFrames] = pupilAnalysis(howManyPlot, displayPlot, createPlot, dFldr, tifFiles, trialNumber, minX, maxX, minY, maxY)
%PUPILANALYSIS Summary of this function goes here
%   Detailed explanation goes here
%     mkdir(dFldr, folderName);
%     disp(['Pupil Analysis Trial ', num2str(trialNumber)])
    %if ~isempty(iTOI)
    %Pupil file info
    %pupilFile = strtrim([dFldr cFNames{iTOI}]);
    fileName = strtrim(tifFiles{trialNumber});
    pupilFile = fullfile(dFldr, fileName);

    pupilInfo = imfinfo(pupilFile);
    nFrames = size(pupilInfo,1);

    
    dataTrialArray = cell(nFrames, 1);
    
    
    % initialize output Figures folder % (can be broken into own file, if needed)
%     if ~exist(fullfile(dFldr, 'Figures'),'dir')
%         mkdir 'Figures'
%     end
    
    tifFileName = sprintf('Trial_%d.tif', trialNumber);
    fullFileName = fullfile(dFldr, 'Figures', tifFileName);
    if exist(fullFileName, 'file')
        delete(fullFileName)
    end
    % end of initialization %
   
    
    for frameNumber = 1:nFrames
        % Load and crop the image from the current frame
        fullImage = imread(pupilFile,frameNumber);
        cropImage=fullImage([minY:maxY], [minX:maxX]);
        
        

        sizeOfData = size(cropImage);
        filteredImage = zeros(sizeOfData);
        for k = 1:sizeOfData(1)
            for j = 1:sizeOfData(2)
                if cropImage(k,j) < 8
                    filteredImage(k,j) = 1;
                end
            end
        end
        skin1 = filteredImage;
        
        skin2 = bwmorph(skin1,'close');
        skin3 = bwmorph(skin2,'open');
        skin4 = bwareaopen(skin3,200);
        skin_test = bwareaopen(skin3,100);
        
%         skin5 = imfill(skin4,'holes');
%         skin6 = bwconvhull(skin5);
%         skin7 = bwconvhull(skin6,'objects');
%         pupilRegion = skin7;
%         [r,c] = find(pupilRegion);
%         CH = convhull(r,c);
        
%         % skin 4
%         p1 = bwconvhull(skin4);
%         o1 = bwconvhull(p1,'objects');
%         [r1,c1] = find(o1);
%         CH1 = convhull(r1,c1);
        
        % skin_test (optimal)
        p2 = bwconvhull(skin_test);
        o2 = bwconvhull(p2,'objects');
        [r2,c2] = find(o2);
        CH2 = convhull(r2,c2);

        % figure(2), hold on
        numPoints = size(CH2, 1);
        oldPoints = zeros(size(CH2,1), 2);
        for i = 1:size(CH2,1)
            oldPoints(i,1) = c2(CH2(i));
            oldPoints(i,2) = r2(CH2(i));
        end
        
        [x_outliers, TF] = rmoutliers(oldPoints(:,1));
        numPoints = numPoints - size(x_outliers, 1);
        points = zeros(numPoints, 2);
        
        oldX = oldPoints(:,1);
        points(:,1) = oldX(TF);
        oldY = oldPoints(:,2);
        points(:,2) = oldY(TF);        
        

        t = (1:numPoints)';
        X = ones(numPoints,3);
        X(:,2) = cos((2*pi)/numPoints*t);
        X(:,3) = sin((2*pi)/numPoints*t);
        y = points(:, 1);
        y = y(:);
        beta = X\y;
        yhat = beta(1)+beta(2)*cos((2*pi)/numPoints*t)+beta(3)*sin((2*pi)/numPoints*t);
%         plot(t,y,'b');
        hold on
%         plot(t,yhat,'r','linewidth',2);
        sin_approximation_x = yhat;


        t = (1:numPoints)';
        X = ones(numPoints,3);
        X(:,2) = cos((2*pi)/numPoints*t);
        X(:,3) = sin((2*pi)/numPoints*t);
        y = points(:, 2);
        y = y(:);
        beta = X\y;
        yhat = beta(1)+beta(2)*cos((2*pi)/numPoints*t)+beta(3)*sin((2*pi)/numPoints*t);
%         plot(t,y,'b');
        hold on
%         plot(t,yhat,'r','linewidth',2);
        
        sin_approximation_y = yhat;
        
        
        
        figure(1), hold on
%         
%         imshow(skin, 'InitialMagnification', 'fit');
%         subplot(3,1,1)
%         imshow(cropImage,  'InitialMagnification', 'fit');
%         colormap parula;
% 
%         axis on;
%         hold on;
%         plot(c(CH), r(CH), "*-", 'Color', 'r');
%        
%         subplot(3,1,2)
%         imshow(skin4,  'InitialMagnification', 'fit');
%         colormap parula;
% 
%         axis on;
%         hold on;
%         plot(c1(CH1), r1(CH1), "*-", 'Color', 'green');
%         
        subplot(2,1,1)
        imshow(cropImage,  'InitialMagnification', 'fit');
        colormap parula;

        axis on;
        hold on;
        plot(c2(CH2), r2(CH2), "*-", 'Color', 'blue');
%         
%         
%         subplot(2,1,2)
%         imshow(skin_test,  'InitialMagnification', 'fit');
%         colormap parula;
% 
%         axis on;
%         hold on;
        plot(sin_approximation_x, sin_approximation_y, "*-", 'Color', 'red');
        hold off
       
        
%         imshow(cropImage,'InitialMagnification','fit');
%         
        
        
        % To get the properties of the pupilRegion call regionprops
        % This returns an object called pupilProperties which has the
        % following fields:
        % 'Centroid': [pupilCenterXCoordinate, pupilCenterYCoordinate]
        % 'MajorAxisLength': length of the major axis of the pupil ellipse
        % 'MinorAxisLength': length of the minor axis of the pupil ellipse
        % 'Area': the number of pixels in the pupilRegion
        % 'BoundingBox': the closest rectangle to the pupilRegion [x, y, width, height] 
%         pupilProperties = regionprops(pupilRegion, 'Centroid','MajorAxisLength','MinorAxisLength', 'Area', 'BoundingBox');
        
        pupilProperties = regionprops(skin_test, 'Centroid','MajorAxisLength','MinorAxisLength', 'Area', 'BoundingBox');
        
        % This piece of code checks whether regionprops was able to
        % identify a pupiilRegion
        N = size(pupilProperties,1);
        if N < 1 || isempty(pupilProperties)
            disp(["pupilProperties was empty for Trial ", trialNumber, " Frame ", frameNumber]);  
            continue
        end
        
        dataTrialArray(frameNumber,1) = {pupilProperties(1)};
        
        % Visualization Protocol
        disp(['Pupil Analysis Trial & Frame ', num2str(trialNumber), ' ', num2str(frameNumber)])
        if createPlot == 1
            if howManyPlot > 0
                displayPlot = "on";
                howManyPlot = howManyPlot - 1;
                disp(["How many plots left ", howManyPlot]);
            else
                displayPlot = "off";
            end
            visualizePupilAnalysis(displayPlot, trialNumber, frameNumber, pupilProperties, fullFileName, fullImage, cropImage, skin, minX, maxX, minY, maxY);
        end
    end
end
