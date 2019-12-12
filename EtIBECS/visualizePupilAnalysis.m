function visualizePupilAnalysis(displayPlot, trialNumber, cnt, pupilProperties, fullFileName, fullImage, cropImage, skin, minX, maxX, minY, maxY)       

% Get Ready for plotting
%     clf;
    
    figure('Visible', displayPlot);
    colormap gray;
    sgtitle(['Trial ' num2str(trialNumber) ', Frame ' num2str(cnt)]);
    
    
% Full Image Plot Top Left Corner
    subplot(221);
    imagesc(fullImage), hold on;
    BoundingBox = pupilProperties(1).BoundingBox;
    rectangle('Position', [minX + BoundingBox(1), minY + BoundingBox(2), BoundingBox(3), BoundingBox(4)], 'EdgeColor', 'r', 'LineWidth', 1);
    title('Full Image');
    
% Cropped Image Plot Bottom Left Corner
    subplot(223);
    axis square
    imagesc(cropImage*10), hold on; %imagesc(fullImage([yStart:yEnd],[xStart:xEnd]));
        % yStart = round(minY + BoundingBox(2));
        % yEnd = round(minY + BoundingBox(2)+ BoundingBox(4));
        % xStart = round(minX + BoundingBox(1));
        % xEnd = round(minX + BoundingBox(1)+ BoundingBox(3));
    title('Cropped Image');
    hold off;
    
% Threshold Image Plot Top Right Corner (Ask April if she wants skin or pupilRegion)
    subplot(222);
    axis square
    imagesc(skin);
    title('Threshold Image');
    
% Pupil Measurements Image Plot Bottom Right Corner
    subplot(224);
    axis square
    imagesc(cropImage*10), hold on;
    
    % Center points
    X = round(pupilProperties(1).Centroid(1), 1);
    Y = round(pupilProperties(1).Centroid(2), 1);
    Area = round(pupilProperties(1).Area, 1);
    % Plotting approximation of pupil based on bounding box and rectangle
    % with full curvature, possibly use ellipse annotation here if can
    % figure out
    rectangle('Position', pupilProperties(1).BoundingBox,'EdgeColor',[1 0 0], 'Curvature', [1,1],'LineWidth',1);
    
    % Plot Center points
    plot(X,Y,'g+');
    text(X - 12,Y - 10,['(',num2str(X),',',num2str(Y),')'],'Color',[0 1 0], 'FontSize', 14);
    text(5,10,['Pupil Area: ',num2str(Area)],'Color',[1 0 0], 'FontSize', 14);
    title('Measurements Image');
    hold off;
    drawnow;
    
    
% Saving to a folder, in a tif stack
    % disp("Started writing figures for trial:");
    % disp(trialNumber);
    % disp("Frame Number:")
    % disp(cnt);
    % 
    
    H = getframe(gcf);
    [X, ~] = frame2im(H);

    imwrite(X, fullFileName,'WriteMode','append');
    

end