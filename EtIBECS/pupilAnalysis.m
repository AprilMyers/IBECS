function [pupilSizes, pupilVids] = pupilAnalysis(dFldr, u, iTOI, cFNames, minX, maxX, minY, maxY)
%PUPILANALYSIS Summary of this function goes here
%   Detailed explanation goes here
    disp(['Pupil Analysis Trial ', num2str(u)])
    if ~isempty(iTOI)
        % Pupil file info
        pupilFile = strtrim([dFldr cFNames{iTOI}]);
        pupilInfo = imfinfo(pupilFile);
        nFrames = size(pupilInfo,1);
        pupilSizeproc = nan(nFrames,1);
        pupilVid = nan(maxY-minY+1,maxX-minX+1,nFrames);
        % Measure pupil
        area_vector = zeros(1,nFrames)
        for cnt = 1:nFrames
            %load and crop
            fullImage = imread(pupilFile,cnt);
            cropImage=fullImage([minY:maxY], [minX:maxX]);
            % Threshold
            skin =~ im2bw(cropImage,0.05);
            skin = bwmorph(skin,'close');
            skin = bwmorph(skin,'open');
            skin = bwareaopen(skin,200);
            skin = imfill(skin,'holes');
            % Measure pupil
            % Select larger area
            L = bwlabel(skin);
            [out_a] = regionprops(L);
            N = size(out_a,1);
            if N < 1 || isempty(out_a) % Returns if no object in the image
                solo_cara=[ ];
                continue
            end
            areas=[out_a.Area];
            [area_max pam]=max(areas);
            % Measure pupil
            centro=round(out_a(pam).Centroid);
            X=centro(1);
            Y=centro(2);
            pupilSizeXY = out_a(pam).BoundingBox;
            sX = pupilSizeXY(3);
            sY = pupilSizeXY(4);
            
            % save data
            pupilSizeproc(cnt) = mean([sX,sY]);
            pupilVid(:,:,cnt) = uint8(cropImage);
            
            % Visualization Stuff
            visualizePupilAnalysis(showAllPlot, showMeasure, u, cnt, fullImage, minX, maxX, minY, maxY, cropImage, skin, pam, X, Y, sX, sY);
        end
        pupilSizes{u} = pupilSizeproc;
        pupilVids{u} = pupilVid;
    end
    
    %Plot pupil size over trial   
    %     figure(2), hold all
    %     plot(pupilSize)
    %     pause
    save([dFldr 'pupilData.mat'],'pupilSizes','pupilVids','-v7.3')
end