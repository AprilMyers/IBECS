function visualizePupilAnalysis(showAllPlot, showMeasure, u, cnt, fullImage, minX, maxX, minY, maxY, cropImage, skin, pam, X, Y, sX, sY)        
    %Visualization Stuff
    if showAllPlot==true | showMeasure==true
    figure(1), clf
    sgtitle(['Trial ' num2str(u) ', Frame ' num2str(cnt)])
    if showAllPlot==true
        % Show full image
        subplot(221),
        imagesc(fullImage), hold on
        xs = [minX maxX maxX minX minX];
        ys = [minY minY maxY maxY minY];
        plot(xs, ys,'r','linestyle','-');
        title('Full')
        % Show crop
        subplot(223)
        imagesc(cropImage*2)
        title('Cropped')
    end
    if showMeasure==true
        % display thresholding
        subplot(222)
        imagesc(skin);
        title('Threshold')
        % Display pupil measurements
        subplot(224)
        title('Measurements')
        imagesc(cropImage*10);
        colormap gray
        hold on
        rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
            'Curvature', [1,1],'LineWidth',1)
        plot(X,Y,'g+')
        text(X+10,Y,['(',num2str(sX),',',num2str(sY),')'],'Color',[1 1 0])
        hold off
    end
    drawnow
    end
end