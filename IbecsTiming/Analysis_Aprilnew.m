



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze IBECS session data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize
%clearvars
%close all

% Session info
dFldrs = {
      'D:\data\M56\20190501_postConditioning\Cond3';...
%          'E:\data\M62\20190511\Cond1';...
%         'E:\data\M62\20190512\Cond2';...
%          'E:\data\M62\20190513\Cond3';...
%         'E:\data\M62\20190515\PostCond2';...
%        'E:\data\M58\20190504\postCond3';...
     
}; % folders containing data
% Other Analsys settings
WL = 1;  %  0/1 = Keep/Remove trials without licks
shkMans = [0 0 0 1 1 0]; % sessions need manual adjustments
shkFlps = {[],[],[],[10 19 28 37 46 55 64 73 82 91],[10 19 28 37 46 55 64 73 82 91],[],[]}; % trials to flip shock on/off
nSesh = length(dFldrs);
% Ryans Figure Settings
RyansFigureSettings
maxLickDiff = [];

%% Lickport analysis
runSeshs = 1:length(dFldrs); % 5; %
for sesh = runSeshs
    dFldr = dFldrs{sesh};
    seshN = 's';  %Session name
    
   
    shkMan = shkMans(sesh); % used a shock protocol but shocker was turned off for some trials
    shkFlp = shkFlps{sesh}; % list of trials to flip shocks on or off
    
    % Files
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
    cFNames = num2cell(ls(dFldr),2); %cell array of file names
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
    seshTitle = cFNames{10}(1:find('_'==cFNames{10},1,'first')-1);

    %% Read in Protocol Settings
    iTOI = fProt;
    fNm = strtrim([dFldr cFNames{iTOI}]);
    protocol = readtable(fNm, 'ReadVariableNames', true, 'Delimiter', '\t');
    
    iTOI = find(fStim);
    fNm = strtrim([dFldr cFNames{iTOI}]);
    stimInfo = load(fNm);
    
    % number of trials
    nTrls = str2num(protocol.TrialName{1});
    nTrlsCheck = length(protocol.TrialName)-1;
    if nTrls~=nTrlsCheck
        warning('Mismatched number of trials in Protocol')
    end
    
    % read trial times (s)
    trlDurs = cellfun(@str2num,protocol.Record(2:end));
    vlvDels = cellfun(@str2num,protocol.FluidValve(2:end));
    vlvDurs = cellfun(@str2num,protocol.Var5(2:end));
    lckDels = cellfun(@str2num,protocol.Lickport(2:end));
    lckDurs = cellfun(@str2num,protocol.Var7(2:end));
    stmDels = stimInfo.staticMovingT(1)*ones(size(trlDurs));
    stmDurs = stimInfo.staticMovingT(2)*ones(size(trlDurs));
    shkDels = cellfun(@str2num,protocol.ShockControl(2:end));
    shkDurs = cellfun(@str2num,protocol.Var10(2:end));        %double check if this is correct
    
    % read from this single trial to get all parameters
    readTrial = 6;
    trlDur = trlDurs(readTrial);
    vlvDel = vlvDels(readTrial);
    vlvDur = vlvDurs(readTrial);
    lckDel = lckDels(readTrial);
    lckDur = lckDurs(readTrial);
    stmDel = stmDels(readTrial);
    stmDur = stmDurs(readTrial);
    
    % Shock times - adjusting for manual changes
    if shkMan==1;
        shkDels(shkFlp) = NaN;
        shkDurs(shkFlp) = NaN;
    end
    shkDels(shkDels==-1) = nan;
    shkDurs(shkDurs==-1) = nan;
    
    shkDel = min(shkDels);
    shkDur = max(shkDurs);
    
    % 0/1 = Don't/Do count spikes that occur during shock time window
    if sum(~isnan(shkDels))>0; shkOn=1; else shkOn=0; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Lick Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Read in Lick data
    
    trls = [1:nTrls];   %trials to process
    lickTimes = cell(size(trls));
    lickCounts = nan(size(trls));
    lickBinCountsT = nan(length(trls),3);
    lickSeq = nan(size(trls))';
    if shkOn==1, shkAdj=shkDur; else shkAdj=0; end %shock adjustment on or off
    
    parfor v = trls %%loop to process trials
        % for v = trls %%loop to process trials
        disp(['Lick Trial_' num2str(v)])
        % select trial
        fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',v)]); %
        
        % Get lick measurement times
        iTOI = find(fSesh & fLick & fSTrl & fTxt);
        if ~isempty(iTOI)
            lickSeq(v) = 1;
            fNm = strtrim([dFldr cFNames{iTOI}]);
            sData = textread(fNm, '%s', 'delimiter', '\n');
            nSamps = length(sData);     % number of samples
            mSecLick = nan(nSamps,1);   % Timestamp of each lick sample
            lLick = nan(nSamps,1);      % Logical of each lick sample
            for u = 1:nSamps
                sStr = sData{u};
                offset = find(sStr=='.',1); % adjusts for change in string length due to more digits
                mSecLick(u) = str2num(sStr(1:offset+6));
                lLick(u) = logical(str2num(sStr(offset+8)));
            end
            sec1Lick = str2num(fNm(end-22:end-15));
            secsLick = mSecLick/1000+ sec1Lick;
            
            % find lick starts
            lickStarts = diff(lLick)==1;
            lickTimes{v} = mSecLick(find(lickStarts)+1)/1000;
            lickCounts(v) = length(lickTimes{v});
        else
            lickSeq(v) = 0;
        end
        % select trial
    end
    
    beep
    
    %% OPTIONAL; Figure out if\where there are changes the sequence
    
    try
        trigTimes = cell(size(nTrls));
        for v = 1:nTrls %%loop to process trials
            %     for v = trls %%loop to process trials
%             disp(['Trial_' num2str(v)])
            % select trial
            fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',v)]); %
            
            % Get lick measurement times
            iTOI = find(fSesh & fTrig & fSTrl & fTxt);
            fNm = strtrim([dFldr cFNames{iTOI}]);
            trigTimes{v} = datetime([fNm(end-26:end-25) ':' fNm(end-23:end-22) ':' fNm(end-20:end-19)]);
        end
        
        for v = 1:nTrls-1
            disp(['Trial_' num2str(v)])
            t1 = trigTimes{v};
            t2 = trigTimes{v+1};
            diffT(v) = seconds(t2-t1);
        end
        disp(diffT')
        unique(diffT)
        
        errorTrials = find(diffT==90) + [1:sum(diffT==90)];
        removeSeq = zeros(1,nTrls);
        removeSeq(errorTrials)= 1;
    end
    %% Trial-by-Trial Licks and PSTH
    
    %initialize
    
    %stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end
    
    %from autolog
    stimSeq0 = zeros(nTrls,1);
    stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
    [~,~,stimId] = unique(stimSeq0);
    stimSeq = stimId-1;
    % stimSeq = zeros(nTrls,1);
    
    % color selection
    nStim = length(unique(stimSeq));
    if nStim==1, colors = {'k','r'};
    else if nStim==2, colors = {'k', 'r'};
        else if nStim==3, colors = {'k', 'b', 'r'};
            else if nStim==5, colors = {'k', 'g', 'b', 'r', 'c'};
                else if nStim==9, colors = {'b','c','g','c','b','m','r','m','k'};
                    end
                end
            end
        end
    end
    colShk = [.8 .9 0];
    colStmLn = 'r';
    
    for v = 1:nTrls % Time-binned lick counts
        
        lickBef =  sum(lickTimes{v}<stmDel);
        lickDur =  sum(lickTimes{v}>stmDel & lickTimes{v}<stmDel+stmDur);
        lickAft =  sum(lickTimes{v}>stmDel+stmDur);
        lickBinCounts(v,:) = [lickBef, lickDur, lickAft];
        
    end
    % remove no-lick trials from grouped analysis?
    
        trialsWLicks = logical(lickBinCounts(:,1));
        lickBinCountsWL = lickBinCounts(trialsWLicks,:);
        stimSeqWL = stimSeq(trialsWLicks);
        
        lickPortTrls = find(lckDels~=-1);
        lickBinCounts0 = lickBinCounts(lickPortTrls,:);
    
    % Lick Rasters and PSTHs
%     figure(1), hold on
    rows = 4;
    cols = nSesh;
    % Licks in Trial Sequence
    subplot(rows,cols,[sesh sesh+nSesh]),  hold all
    title(seshTitle)
%     xlabel('Time (s)')
 if sesh==1, ylabel('Trial'); end
    line([0 trlDur],[0 0],'color','k')
    line([stmDel stmDel],[0 v+1],'color',colStmLn)
    line([stmDel+stmDur stmDel+stmDur],[0 v+1],'color',colStmLn)
    for v = 1:nTrls
        %draw shock time
        if shkOn==1; trlShkDel = shkDels(v);
            if trlShkDel>=0; trlShkDur = shkDurs(v);
                line([trlShkDel trlShkDel+trlShkDur],[v v],'color',colShk,'linewidth',6)
            end
        end
        
        %draw lick times
        lickTime = lickTimes{v};
        color = colors{stimSeq(v)+1};
        scatter(lickTime, ones(size(lickTime))*v,color,'.')
    end
    axis tight
    
    % Calculate PSTHs
    binS = .25;
    edges = 0:binS:trlDur;
    xTime = edges(1:end-1)+binS/2;
    psthTrials = nan(nTrls,length(edges)-1);
    for v = 1:nTrls
        lickTime = lickTimes{v};
        if isempty(lickTime), lickTime=zeros(0,1); end
        psthTrial = histc(lickTime,edges)/binS;
        psthTrials(v,:) = psthTrial(1:end-1);
    end
    windLength = 5;
    smthFlt = gausswin(windLength);
    smthFltS = smthFlt/sum(smthFlt);
    % Plot PSTH
    subplot(rows,cols,2*nSesh+sesh), hold on
%     title({'PSTH Lick Rate, Lickless Trials Removed'})
    xlabel('Time (s)')
   if sesh==1, ylabel('Lick rate (Hz)'); end
    maxY = 0;
    for v = 1:nStim
        s = v-1;
        mS = conv(mean(psthTrials(stimSeq==s&lickSeq&trialsWLicks,:)),smthFltS,'same');
        eS = conv(ste(psthTrials(stimSeq==s&lickSeq&trialsWLicks,:)),smthFltS,'same');
        color = colors{v};
        patch([xTime fliplr(xTime)], [mS+eS fliplr(mS-eS)],color,'FaceAlpha',.1,'edgecolor','none')
        plot(xTime, mS,color,'linewidth',2)
        maxY = max([maxY mS+eS]);
    end
    if shkOn==1;
        line(repmat(min(shkDels),1,2),[0 maxY],'color',colShk,'linewidth',2);
        line(repmat(sum(max([shkDels shkDurs])),1,2),[0 maxY],'color',colShk,'linewidth',2);
    end
    line([stmDel stmDel],[0 maxY],'linestyle','--','linewidth',2,'color',colStmLn)
    line([stmDel+stmDur stmDel+stmDur],[0 maxY],'linestyle','--','linewidth',2,'color',colStmLn)
    axis tight
    
    %% blank:stim ratio
    stimTrls = find(stimSeq);
    stimLickTrls = stimTrls(1:2:end)/2;
    befTrls = stimLickTrls(1:end,:)-1;
    aftTrls = stimLickTrls(1:end,:)+1; %blank trials flanking stimulus trials
    flankLicks = mean([lickBinCounts0(befTrls,2) lickBinCounts0(aftTrls,2)] ,2);
    stimLicks = lickBinCounts0(stimLickTrls(1:end),2);
    stimFlankDif = stimLicks-flankLicks;
    noLickTrls = floor((find(lickBinCounts0(:,1)==0)+2)/3);             % find no-lick trials
    shkTrls = round(find(lckDels==2&stimSeq&~isnan(shkDels))/2/3);  % find shock trials
    stimFlankDif([noLickTrls; shkTrls]) = NaN;                      % remove no-lick trials and shock trials
    maxLickDiff = max([maxLickDiff max(abs(stimFlankDif))]);
    subplot(rows,cols,nSesh*3+sesh)
    bar(stimFlankDif,'facecolor','r'), box off
    if sesh==1, ylabel('Lick Count Difference'); end
    xlabel('Stim & Lick Trials (in sequence)')
    [H P(sesh)] = ttest(stimFlankDif,0,'tail','left');
    

    
    %%
end
for sesh = 1:nSesh
subplot(rows,cols,nSesh*3+sesh)
ylim([-maxLickDiff maxLickDiff])
text(1,maxLickDiff/2,['p = ' num2str(P(sesh))],'FontSize', 18)
end
%% group differences
figure(2), clf
rows = 1;
cols = 4;

lickRate = lickBinCountsWL;
lickRate(:,1) = lickRate(:,1)/(stmDel-lckDel-.35);
lickRate(:,2) = lickRate(:,2)/stmDur;
lickRate(:,3) = lickRate(:,3)/(lckDur+lckDel-stmDel-stmDur);
lickRatio = lickRate(:,1)./lickRate(:,2);
lickRatio(lickRatio==Inf)=nan;

% lick rate before, during after stimulus
subplot(rows,cols,1),  hold all
title('lick rate before, during after stimulus')
boxplot(lickRate)

set(gca,'xtick',[1:3],'xticklabel',{'Before','During','After'})
ylabel('Lick Rate')
figure, [~,~,STATS] = anova1(lickRate);
multcompare(STATS)

% lick rate during stimulus vs blank
figure(2), hold all, subplot(rows,cols,2), hold all
title('lick rate during stimulus vs blank')
boxplot(lickRate(:,2), stimSeqWL>0)
set(gca,'xtick',[1:2],'xticklabel',{'Blank','Stimulus'})
figure,[~,~,STATS] = anova1(lickRate(:,2), stimSeqWL>0,'off');
[COMPARISON,MEANS,H,GNAMES]=multcompare(STATS);
COMPARISON(end)

% lick rate during stimulus
figure(2), hold all,  subplot(rows,cols,3), hold all
title('lick rate during stimulus')
boxplot(lickRate(:,2), stimSeqWL)
set(gca,'xtick',[1:5],'xticklabel',{'Blank','Down','Left','Up','Right'})
figure, [~,~,STATS] = anova1(lickRate(:,2), stimSeqWL,'off');
[COMPARISON,MEANS,H,GNAMES]=multcompare(STATS);
multcompare(STATS)

% lick ratio before:during stimulus
figure(2), hold all,  subplot(rows,cols,4), hold all
title('lick ratio before:during stimulus')
boxplot(lickRatio, stimSeqWL)
set(gca,'xtick',[1:5],'xticklabel',{'Blank','Down','Left','Up','Right'})
ylabel('Lick Ratio (Before/During')
figure, [~,~,STATS] = anova1(lickRatio, stimSeqWL);
[COMPARISON,MEANS,H,GNAMES]=multcompare(STATS);
COMPARISON(end)
multcompare(STATS)
% [h,p]=ttest2(lickRate(stimSeq==0,2), lickRate(stimSeq==3,2))
% [h,p]=ttest2(lickRate(stimSeq==0,2), lickRate(stimSeq==4,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pupil Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Eye

% select trial 1
fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',1)]);
% Get first cam1 trial video
iTOI = find(fSesh & fCam1 & fTiff & fSTrl);
pupilFile = strtrim([dFldr cFNames{iTOI}]);
           
fullImage = imread(pupilFile,1);
figure

acceptCrop=0
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
acceptCrop= input('Accept Crop? Enter 1 to accept or 0 to reject: ')
end
%% Measure Pupil Session

%plot settings

showAllPlot = false;
showMeasure = false;

startFromTrial = 1;
figure(2), clf
tic
parfor u = startFromTrial:nTrls
   %pause
    disp(['Pupil Analysis Trial ', num2str(u)])
    % select trial
    fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',u)]); %
    % Get pupil measurements
    iTOI = find(fSesh & fCam1 & fSTrl & fTiff);
    
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
            cropImage=fullImage([minY:maxY],[minX:maxX]);
            % Threshold
            skin=~im2bw(cropImage,0.05);
            skin=bwmorph(skin,'close');
            skin=bwmorph(skin,'open');
            skin=bwareaopen(skin,200);
            skin=imfill(skin,'holes');
            % Measure pupil
            % Select larger area
            L=bwlabel(skin);
            [out_a]=regionprops(L);
            N=size(out_a,1);
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
        pupilSizes{u} = pupilSizeproc;
        pupilVids{u} = pupilVid;
    end
    
    %Plot pupil size over trial   
    %     figure(2), hold all
    %     plot(pupilSize)
    %     pause
end
toc
beep

save([dFldr 'pupilData.mat'],'pupilSizes','pupilVids','-v7.3')



%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
for u = 1:nTrls
    pupilSize = pupilSizes{u};
    
    pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
end

pupilAbsMax = max(pupilSize500(:));
% pupilAbsMax = max(abs(pupilSize500(:)));
pupilAllNorm = pupilSize500/pupilAbsMax;
% wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
    % compare raw, interpolated and smoothed data
    % subplot(rows,cols,1), cla, hold on
    % line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
    % plot(wheelTimes{v},wheelCounts{v},'-','color',color)
    % plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
    % plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
    % axis tight, xlim([0 trlDur]),  box off
    % pause
    color = colors{trialType(v)};
    pupilTimeFrame = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
    plot(pupilTimeFrame,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,pupilAllNorm(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,pupilAllNorm(typeTrials,:),'-','color',color)
    ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(pupilAllNorm(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((pupilAllNorm(typeTrials,:)));
    pupilSte = nanste((pupilAllNorm(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,pupilAllNorm(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
pupilSizes{u} = pupilSize
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Loop through to fill emptycell with pupilareaProc data

pupilFile = strtrim([dFldr cFNames{iTOI}]);
pupilInfo = imfinfo(pupilFile);
nFrames = size(pupilInfo,1);
pupilProc = proc.pupil.area;
pupilareaProc = cell(1,nTrls);

framebegin = 1;
frameend = nFrames;

for u = startFromTrial:nTrls
    fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',u)]); %
    iTOI = find(fSesh & fCam1 & fSTrl & fTiff);
    pupilFile = strtrim([dFldr cFNames{iTOI}]);
    pupilInfo = imfinfo(pupilFile);
    nFrames = size(pupilInfo,1);
    frameend = framebegin + nFrames - 1;
    pupilareaProc{u} = pupilProc(framebegin:frameend);
    framebegin = framebegin + nFrames;
end

%// pad columns with necessary number of NaN
b = cellfun( @(c) [c(:) ; NaN(501-numel(c),1)], pupilareaProc ,'un',0 ) ;
%// reshape to get final desired result
pupilareaProc = b(:).';
%pupilareaProc = smoothdata(cell2mat(pupilareaProc), 'gaussian');
%pupilareaProc = mat2cell(pupilareaProc, 501, [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
pupilareaProc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proc pupil plots
figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog 
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
for u = 1:nTrls
    pupilSizeproc = pupilareaProc{u};
    
    pupilSize500proc(u,1:length(pupilSizeproc)) = (pupilSizeproc-pupilSizeproc(1))./pupilSizeproc;
end

pupilAbsMaxproc = max(pupilSize500proc(:));
% pupilAbsMax = max(abs(pupilSize500(:)));
pupilAllNormproc = pupilSize500proc/pupilAbsMaxproc;
% wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
    % compare raw, interpolated and smoothed data
    % subplot(rows,cols,1), cla, hold on
    % line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
    % plot(wheelTimes{v},wheelCounts{v},'-','color',color)
    % plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
    % plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
    % axis tight, xlim([0 trlDur]),  box off
    % pause
    color = colors{trialType(v)};
    pupilTimeFrame = 1:length(pupilSize500proc(v,:));

    subplot(5,4,4*4+[1 2]), hold all
    plot(pupilTimeFrame,pupilSize500proc(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,pupilAllNormproc(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,pupilAllNormproc(typeTrials,:),'-','color',color)
    ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(pupilAllNormproc(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMeanproc = nanmean((pupilAllNormproc(typeTrials,:)));
    pupilSteproc = nanste((pupilAllNormproc(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMeanproc+pupilSteproc fliplr(pupilMeanproc-pupilSteproc)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMeanproc,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,pupilAllNormproc(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
pupilareaProc{u} = pupilSizeproc
hold off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional pupil plots

xs= [1:length(lickCounts)];
%plt 1
rows = 3;
cols = 1;
figure(3), clf, hold on

subplot(rows,cols,1), hold all
title('Trial total lick count in seuqnce')
ylabel('Lick Count')
plot(lickCounts)
scatter(xs,lickCounts)
lsline

%plt 2
subplot(rows,cols,2), hold all
title('Binned lick rate in seuqnce')
ylabel('lick rate per bin')
plot(lickRate)
legend('Before','During','After')
% trialXs = repmat(xs,3,1)';
% scatter(trialXs(:), lickRate(:))
% lsline

%plt 3
subplot(rows,cols,3), hold all
title('Binned lick ratio, Before:During, in seuqnce')
xlabel('Time(s)')
ylabel('Ratio of licks, Before:DuringLick')
plot(lickRatio)
scatter(1:length(lickRatio),lickRatio(:))
lsline

%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proc Additional pupil plots




pupilMeanarea_ses = cell(1,91)
xs= [1:length(pupilMeanarea_ses)];
for u = startFromTrial:nTrls ; 
    pupilproc_u = pupilareaProc{u} ;
    pupilMeanarea_ses{u} = mean(pupilproc_u)
end
   
%plt 1
rows = 2;
cols = 1;
figure(2), clf, hold on

subplot(rows,cols,1), hold all
title('Mean pupil area for each session')
ylabel('Pupil area')
plot(cell2mat(pupilMeanarea_ses))
scatter(xs,(cell2mat(pupilMeanarea_ses)))
lsline

% %plt 2    
subplot(rows,cols,2), hold all
title('Pupil areas plotted for all session')
ylabel('Pupil area')
plot(proc.pupil.area)
% scatter(xs,(cell2mat(pupilMeanarea_ses)))
lsline

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean data
% Compute list of instantaneous dt for all sessions
% Filter pupilareaProc based on pupilareadt

pupilareadt = cell(1,91) ;
for u = startFromTrial:nTrls ; 
    pupilproc_u = pupilareaProc{u} ;
    t = 1:length(pupilproc_u) ;
    v = zeros(length(t)-1,1) ;
    for i = 1:length(t)-1 ;
      v(i) = (pupilproc_u(i+1)-pupilproc_u(i))/(t(i+1)-t(i)) ;
      pupilareadt{u} = v ;
    end
end
pupilareadt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % Filter pupilareaProc based on pupilareadt
  
% area_ses1_dt = single(pupilareadt{1}) ;
% area_ses1_dt = 
% 
% ntrials = 91
% length = 500
%     for trialindex
%         areases1(:,index) = pupilareareadt{1,trialindex}(:)
%     end
% end

area_ses_done = cell(1,91);
for index=1:nTrls
    disp(index)
    area_ses1_dt = single(pupilareadt{index}) ;
    area_ses1_dt = [NaN;pupilareadt{index}(1:end)];
    
    area_ses1 = pupilareaProc{index};
    %ses1 = table(area_ses1_dt, area_ses1);
    % extract first row: ses1(1:501, 1)
    isgoodframe = (-50 < area_ses1_dt & area_ses1_dt < 50) ;
    area_ses1(~isgoodframe) = nan
    
    area_ses_done{index} = area_ses1
end

pupilareaProc = area_ses_done;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize pupil session 

startFromTrial = 1;
areaMatrix = zeros(501, nTrls);
for u = startFromTrial:nTrls
    % Pupil file info
    pupilVid=pupilVids{u};
    nFrames = length(pupilVid);
    
    areaVec = zeros(501,1)
    
    % Measure pupil
    for cnt = 1:nFrames
        disp(num2str(cnt))
        
        % load and crop
        cropImage=uint8(pupilVid(:,:,cnt));
        % Threshold
        skin=~im2bw(cropImage,0.05);
        skin=bwmorph(skin,'close');
        skin=bwmorph(skin,'open');
        skin=bwareaopen(skin,200);
        skin=imfill(skin,'holes');
        % Measure   
        % Select larger area
        L=bwlabel(skin);
        
        out_a = regionprops(L);
        area =regionprops(L,'Area');
    
        areaVec(cnt) = area.Area;

        N=size(out_a,1);
        if N < 1 || isempty(out_a) % Returns if no object in the image
            solo_cara=[ ];
            continue
        end
            %areas=[out_a.Area];
            %[area_max pam]=max(areas);
        % Measure pupil
        centro=round(out_a(pam).Centroid);
        X=centro(1);
        Y=centro(2);
        pupilSizeXY = out_a(pam).BoundingBox;
        sX = pupilSizeXY(3);
        sY = pupilSizeXY(4);
        
        
        figure(1), clf
        sgtitle(['Trial ' num2str(u) ', Frame ' num2str(cnt)])
        
        % display thresholding
        %         subplot(121)
        %         imagesc(skin);
        %         title('Threshold')
        % Display pupil measurements
        subplot(122)
        title('Measurements')
        imagesc(cropImage*10);
        colormap gray
        hold on
        rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
            'Curvature', [1,1],'LineWidth',1)
        plot(X,Y,'g+')
        text(X+10,Y,['(',num2str(sX),',',num2str(sY),')'],'Color',[1 1 0])
        hold off
        
%         drawnow
    end
    areaMatrix(:,u) = areaVec;
    %    pause
end
toc
beep   
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize pupil session (filter w/ Proc)
 
%getPupilVideoFromProc(proc)
startFromTrial = 1;
xs = cell(nTrls,1);
ys = cell(nTrls,1);
for u = startFromTrial:nTrls
    % Pupil file info
    pupilVid=pupilVids{u};
    nFrames = length(pupilVid);
    file = proc.files{u};
    disp("The file is:")
    disp(file)
    x0 = round(proc.locROI{1,5}(1) * proc.sc) ;
    y0 = round(proc.locROI{1,5}(2) * proc.sc);
    x1 = round(x0 + proc.locROI{1,5}(3) * proc.sc);
    y1 = round(y0 +  proc.locROI{1,5}(4) * proc.sc);
    disp("x0: " + x0)
    disp("y0: " + y0)
    disp("x1: " + x1)
    disp("y1: " + y1)
    
    mov = VideoReader(file);
    nFrames = mov.NumberOfFrames;
    %end ethan code
    
    
    % Measure pupil
    for cnt = 1:nFrames

        disp(num2str(cnt))
         % load and crop
        cropImage=uint8(pupilVid(:,:,cnt));
        % Threshold
        skin=~im2bw(cropImage,0.05);
        skin=bwmorph(skin,'close');
        skin=bwmorph(skin,'open');
        skin=bwareaopen(skin,200);
        skin=imfill(skin,'holes');
        % Measure
        
        % Select larger area
        L=bwlabel(skin);
        out_a=regionprops(L);
        N=size(out_a,1);
        if N < 1 || isempty(out_a) % Returns if no object in the image
            solo_cara=[ ];
            continue
        end   
        areas=[out_a.Area];
        [area_max pam]=max(areas);
        % Measure pupil
        centro=out_a(pam).Centroid;%round(out_a(pam).Centroid);
        X=centro(1)-1;
        Y=centro(2)-1;
        pupilSizeXY = out_a(pam).BoundingBox
        sX = pupilSizeXY(3);
        sY = pupilSizeXY(4);
        
        % start ethan code for bringning in proc pupil data to plot ellipse
        % remember to account for the nframes later
        centerOfPupilX = proc.pupil.com(cnt,1);
        centerOfPupilY = proc.pupil.com(cnt,2); 
        
        
        figure(1), clf
        sgtitle(['Trial ' num2str(u) ', Frame ' num2str(cnt)])
        
%         display thresholding
%                 subplot(121)
%                 imagesc(skin);
%                 title('Threshold')
%         Display pupil measurements
%         subplot(122)
        title('Measurements')
        imagesc(cropImage*10);
        colormap gray
        hold on
        temp = pupilareaProc(1,  int8(u));
        temp2 = cellfun(@isnan,temp, 'UniformOutput', false)
        
        if temp2{1,1}(cnt) == 0
            rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],'Curvature', [1,1],'LineWidth',1)
        end
        %ethan plotting center of pupil from proc
        %plot(centerOfPupilX, centerOfPupilY,'b+')
        % end ethan plotting
        plot(X,Y,'g+')
        text(X+10,Y,['(',num2str(sX),',',num2str(sY),')'],'Color',[1 1 0])
        hold off
        
        %save X and Y
        xs{u}(cnt) = X;
        ys{u}(cnt) = Y; 
        pause
    end
end
toc
beep   


%% troubleshooting
figure, hold all
badTrials = [13];
for u = badTrials(end):nTrls
    disp(num2str(u))
    plot(pupilSizes{u})
    pause
end