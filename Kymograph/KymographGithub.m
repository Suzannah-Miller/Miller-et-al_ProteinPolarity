%% Load movie
clear 
close all

% Path for bioformats toolbox
addpath '/Applications/MATLAB_R2020a.app/toolbox/bfmatlab' 

% Directory for .tif movies
dir = '/Volumes/Seagate4Tb/MovieAnalysis/NikonWidefield/CroppedMovies/Moesin/';
cd(dir)

% Name of movie
imname = '20220502_fitc-180ms_txred_180ms_dapi_500ms_40x_3well001_xy11-Crop';
imgbf = bfopen([imname, '.tif']);
mov = imgbf{1,1};

idxCh1 = find(contains(string(mov(:,2)), "C?=1/3")); % LifeAct-mTagBFP
idxCh2 = find(contains(string(mov(:,2)), "C?=2/3")); % MCAM-GFP
idxCh3 = find(contains(string(mov(:,2)), "C?=3/3")); % ERM-mCherry

% Load movie
for f=1:length(idxCh2)
    movCh1(:,:,f) = im2double(mov{idxCh1(f),1});
    movCh2(:,:,f) = im2double(mov{idxCh2(f),1});
    movCh3(:,:,f) = im2double(mov{idxCh3(f),1});
end
%% Make kymograph and background-subtracted movie
% This code shows the contrast and colormap used to generate Figure 8.
% For quantification of polarization phenotypes, contrast was linearly scaled
% between the 0.1th and 99th percentile intensities, and kymographs were
% saved using the "jet" colormap.

close all
clearvars -except movCh1 movCh2 movCh3 imname dir

cd(dir)

% Name of binary mask of the line along which to plot kymograph (made in
% ImageJ)
linename = [imname, '_Line_5'];
lineMask = imbinarize(imread([linename, '.tif']));
lineStats = regionprops(lineMask, 'Orientation', 'PixelIdxList');

% Angle to rotate cell to the horizontal
angle = -lineStats.Orientation;
imshow(imrotate(lineMask, angle)) % Make sure line is on the right side of image. Add 180 to angle if it is horizontal on the left side.

% Make a mask of the line rotated to the horizontal which is 10 pixels in
% height
[row, col] = find(imrotate(lineMask, angle) == 1);
colStart = min(col); % Location on the rotated image where the kymograph will start
rowStart = round(mean(row));
lineRotated = zeros(size(imrotate(lineMask, angle)));
lineRotated(((rowStart-4):(rowStart+5)),((colStart:size(lineRotated,2)))) = 1; % Thicken line

% Indicates what movie frames to plot
frames = 11:91;

% Thresholds for blobMask
thresh = 0.0057; % Threshold for ch2 (MCAM) for masking - adjust on per image basis
thresh3 = 0.003; % Threshold for ch3 (ERM) for masking - adjust on per image basis

% Make each row of kymograph (find intensities for each frame of movie).
for f = frames
    % Load a single movie frame
    frameCh1 = movCh1(:,:,f);
    frameCh2 = movCh2(:,:,f);
    frameCh3 = movCh3(:,:,f);
    
    % Create mask of the cells in the frame for channel 2 (MCAM-GFP)
    [mask2, blob, imgFilt] = blobMask(frameCh2, thresh);
    statsMask2 = regionprops(mask2, 'Area', 'PixelIdxList');
    maskMain2 = zeros(size(mask2));
    maskMain2(statsMask2(find([statsMask2.Area] == max([statsMask2.Area]))).PixelIdxList) = 1;
    
    % Create mask of the cells in the frame for channel 3 (ERM-mCherry)
    [mask3, blob, imgFilt] = blobMask(frameCh3, thresh3);
    statsMask3 = regionprops(mask3, 'Area', 'PixelIdxList');
    maskMain3 = zeros(size(mask3));
    maskMain3(statsMask3(find([statsMask3.Area] == max([statsMask3.Area]))).PixelIdxList) = 1;
    
    % Find the location of the edge of the original image along the
    % kymograph line
    cropMask = imrotate(ones(size(frameCh1)), angle) + 2.*lineRotated;
    cropIdx = find(cropMask == 2);
    [rowCrop, colCrop] = ind2sub(size(cropMask), cropIdx);
    colEnd = min(colCrop);

    % Mask for the cell of interest
    maskMain = double((maskMain2 + maskMain3) > 0);
        
    % Mask of the region of the cell measured by kymograph 
    % (based on the cell mask). Used to find contrast limits.
    maskRotated = imrotate(maskMain, angle);
    maskRotated(lineRotated == 0) = NaN;
    kymoMaskTemp = mean(maskRotated(:, colStart:colEnd), 1, 'omitnan');
    kymoMaskFilled(f,:) = kymoMaskTemp;
    
    % Combine masks from both channel and dilate to make mask for
    % background subtraction
    maskBg = imdilate((mask3 + mask2) > 0, strel('disk', 3));
    
    % Get values for background subtraction
    bgCh1 = median(frameCh1(maskBg == 0), 'all');
    bgCh2 = median(frameCh2(maskBg == 0), 'all');
    bgCh3 = median(frameCh3(maskBg == 0), 'all');
    
    % Collect background subtraction values and mean signal for each frame
    % (for quality control)
    bg(f,:) = [bgCh1 bgCh2 bgCh3];
    sig(f,:) = [mean(frameCh1(maskMain == 1),'all') mean(frameCh2(maskMain == 1),'all') mean(frameCh3(maskMain == 1),'all')];
    
    % Background-subtracted image
    frameCh1Bg = frameCh1 - bgCh1;
    frameCh1Bg(frameCh1Bg < 0) = 0;
    frameCh2Bg = frameCh2 - bgCh2;
    frameCh2Bg(frameCh2Bg < 0) = 0;
    frameCh3Bg = frameCh3 - bgCh3;
    frameCh3Bg(frameCh3Bg < 0) = 0;
    
    % Rotate images and make kymograph
    frameCh1Rotated = imrotate(frameCh1Bg, angle);
    frameCh1Rotated(lineRotated == 0) = NaN;
    kymoCh1(f,:) = mean(frameCh1Rotated(:, colStart:colEnd), 1, 'omitnan');
    
    frameCh2Rotated = imrotate(frameCh2Bg, angle);
    frameCh2Rotated(lineRotated == 0) = NaN;
    kymoCh2(f,:) = mean(frameCh2Rotated(:, colStart:colEnd), 1, 'omitnan');
    
    frameCh3Rotated = imrotate(frameCh3Bg, angle);
    frameCh3Rotated(lineRotated == 0) = NaN;
    kymoCh3(f,:) = mean(frameCh3Rotated(:, colStart:colEnd), 1, 'omitnan');
    
    % Add frame to background-subtracted movie
    movBgCh1(:,:,f) = frameCh1Bg;
    movBgCh2(:,:,f) = frameCh2Bg;
    movBgCh3(:,:,f) = frameCh3Bg;
end

% Remove frames at the beginning of movie which were omitted (if any)
bg(1:frames(1)-1,:) = [];
sig(1:frames(1)-1,:) = [];
kymoCh1(1:frames(1)-1,:) = [];
kymoCh2(1:frames(1)-1,:) = [];
kymoCh3(1:frames(1)-1,:) = [];
kymoMaskFilled(1:frames(1)-1,:) = [];

% Plot signal, background-subtracted signal, and median background signal
% for each frame to make sure background-subtraction is not doing something
% weird.
figure
subplot(1,3,1) % Median background intensity for each channel
plot(1:length(frames), bg(:,1), 'blue') % Channel 1
hold on
plot(1:length(frames), bg(:,2), 'green') % Channel 2
plot(1:length(frames), bg(:,3), 'red') % Channel 3
hold off
subplot(1,3,2) % Background-subtracted signal
plot(1:length(frames), sig(:,1)-bg(:,1), 'blue')
hold on
plot(1:length(frames), sig(:,2)-bg(:,2), 'green')
plot(1:length(frames), sig(:,3)-bg(:,3), 'red')
hold off
subplot(1,3,3) % Raw signal
plot(1:length(frames), sig(:,1), 'blue')
hold on
plot(1:length(frames), sig(:,2), 'green')
plot(1:length(frames), sig(:,3), 'red')
hold off

% Calculate the upper and lower limits to linearly scale the contrast
scaleCh1 = [prctile(kymoCh1(kymoMaskFilled == 1), .2) 1.4*prctile(kymoCh1(kymoMaskFilled == 1),99)];
scaleCh2 = [prctile(kymoCh2(kymoMaskFilled == 1), .07) 1.4*prctile(kymoCh2(kymoMaskFilled == 1),99)];
scaleCh3 = [min(kymoCh3(kymoMaskFilled == 1)) 1.4*prctile(kymoCh3(kymoMaskFilled == 1),99)];

% Original scale used for evaluation of polarization to the cell periphery.
% scaleCh1 = [prctile(kymoCh1(kymoMaskFilled == 1), .1) prctile(kymoCh1(kymoMaskFilled == 1),99)];
% scaleCh2 = [prctile(kymoCh2(kymoMaskFilled == 1), .1) prctile(kymoCh2(kymoMaskFilled == 1),99)];
% scaleCh3 = [prctile(kymoCh3(kymoMaskFilled == 1), .1) prctile(kymoCh3(kymoMaskFilled == 1),99)];

% Number how long in pixels each frame should be on the kymograph
tWidth = 20;

% Adjust the contrast of the kymograph and movie frames for visualization
for t = frames
    f = t-frames(1)+1;

    kymoShowCh1((f-1)*tWidth+1:(f-1)*tWidth+tWidth,:) = repmat((kymoCh1(f,:)-scaleCh1(1))/(scaleCh1(2)-scaleCh1(1)),tWidth,1);
    kymoShowCh2((f-1)*tWidth+1:(f-1)*tWidth+tWidth,:) = repmat((kymoCh2(f,:)-scaleCh2(1))/(scaleCh2(2)-scaleCh2(1)),tWidth,1);
    kymoShowCh3((f-1)*tWidth+1:(f-1)*tWidth+tWidth,:) = repmat((kymoCh3(f,:)-scaleCh3(1))/(scaleCh3(2)-scaleCh3(1)),tWidth,1);
    
    moveTempCh1 = (movBgCh1(:,:,t) - scaleCh1(1))./(scaleCh1(2)-scaleCh1(1));
    moveTempCh2 = (movBgCh2(:,:,t) - scaleCh2(1))./(scaleCh2(2)-scaleCh2(1));
    moveTempCh3 = (movBgCh3(:,:,t) - scaleCh3(1))./(scaleCh3(2)-scaleCh3(1));
    movBgScaledCh1(:,:,t) = moveTempCh1;
    movBgScaledCh2(:,:,t) = moveTempCh2;
    movBgScaledCh3(:,:,t) = moveTempCh3;
    
    % Make a movie with the kymograph line overlaid
    moveTempCh1(lineMask) = 1;
    moveTempCh2(lineMask) = 1;
    moveTempCh3(lineMask) = 1;
    movBgScaledCh1Line(:,:,t) = moveTempCh1;
    movBgScaledCh2Line(:,:,t) = moveTempCh2;
    movBgScaledCh3Line(:,:,t) = moveTempCh3;
end

% Figure showing the location of kymograph line on the first frame of the
% movie before and after rotation
figure
subplot(1,2,1)
imshowpair(imadjust(movCh2(:,:,frames(1)), scaleCh2+bgCh2), imdilate(lineMask, strel('disk',5)))
subplot(1,2,2)
imshowpair(imrotate(imadjust(movCh2(:,:,frames(1)), scaleCh2+bgCh2), angle), lineRotated)

% Save kymographs and background-subtracted movies
cd('/Users/suzannah/Library/CloudStorage/OneDrive-UCB-O365/AhnLabResearch/Experiments/MovieAnalysis/Kymographs/Figures/Paper/')

% % False-colored with pink colormap
map = colormap(pink(256)); % For quantification of polarization phenotypes, colormap(jet(256)) was used.
% imwrite(ind2rgb(uint8(kymoShowCh1*255),map), [linename, '_kymoCh1_pink_lowerC_Long.tif'], 'Compression', 'none')
% imwrite(ind2rgb(uint8(kymoShowCh2*255),map), [linename, '_kymoCh2_pink_lowerC_Long.tif'], 'Compression', 'none')
% imwrite(ind2rgb(uint8(kymoShowCh3*255),map), [linename, '_kymoCh3_pink_lowerC_Long.tif'], 'Compression', 'none')

% % Grayscale
% imwrite(uint8(kymoShowCh1*255), [linename, '_kymoCh1_grey_lowerC.tif'], 'Compression', 'none')
% imwrite(uint8(kymoShowCh2*255), [linename, '_kymoCh2_grey_lowerC.tif'], 'Compression', 'none')
% imwrite(uint8(kymoShowCh3*255), [linename, '_kymoCh3_grey_lowerC.tif'], 'Compression', 'none')

cd('/Users/suzannah/Library/CloudStorage/OneDrive-UCB-O365/AhnLabResearch/Experiments/MovieAnalysis/Kymographs/Figures/Paper/')
map = colormap(pink(256));
for f = frames
    % Background-subtracted movie (false colored)
%     imwrite(ind2rgb(uint8(movBgScaledCh1(:,:,f).*255),map), [imname, 'ScaledMovieCh1_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%      imwrite(ind2rgb(uint8(movBgScaledCh2(:,:,f).*255),map), [imname, 'ScaledMovieCh2_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%     imwrite(ind2rgb(uint8(movBgScaledCh3(:,:,f).*255),map), [imname, 'ScaledMovieCh3_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%   
    % Background-subtracted movie (false colored) with kymograph line
    % overlaid
    %imwrite(ind2rgb(uint8(movBgScaledCh1Line(:,:,f).*255),map), [linename, 'ScaledMovieCh1_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%      imwrite(ind2rgb(uint8(movBgScaledCh2Line(:,:,f).*255),map), [linename, 'ScaledMovieCh2_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%     imwrite(ind2rgb(uint8(movBgScaledCh3Line(:,:,f).*255),map), [linename, 'ScaledMovieCh3_pink_lowerC.tif'], 'Compression', 'none','WriteMode','append')
end

% % Background-subtracted movie (greyscale)
% for f = frames
%     imwrite(uint8(movBgScaledCh1(:,:,f).*255), [imname, '_ScaledMovieCh1_Grey_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%     imwrite(uint8(movBgScaledCh2(:,:,f).*255), [imname, '_ScaledMovieCh2_Grey_lowerC.tif'], 'Compression', 'none','WriteMode','append')
%     imwrite(uint8(movBgScaledCh3(:,:,f).*255), [imname, '_ScaledMovieCh3_Grey_lowerC.tif'], 'Compression', 'none','WriteMode','append')
% end

% Check thresholding of cell (used for background subtraction and marking
% edge of cell)
cd(dir)
bFrames = struct('B', {});
for f = frames([5 20 77])
    frameCh2 = movCh2(:,:,f);
    [mask, blob, imgFilt] = blobMask(frameCh2, thresh);
    statsMask = regionprops(mask, 'Area', 'PixelIdxList');
    maskMain = zeros(size(mask));
    maskMain(statsMask(find([statsMask.Area] == max([statsMask.Area]))).PixelIdxList) = 1;
    B = bwboundaries(maskMain);
    maskFrames(:,:,f) = maskMain;
    bFrames(f).B = B{1};
    
    pctl = prctile(frameCh2(maskMain == 1), 70);
    pctlMask = (frameCh2 > pctl).*maskMain;
    
    figure
    subplot(1,2,1)
    imshowpair(imadjust(movCh2(:,:,f), scaleCh2), maskMain)
    hold on
    boundary = B{1};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', .5)
    hold off
    subplot(1,2,2)
    imshow(maskMain+pctlMask,[0 2])
    
end

cd(dir)
bFrames = struct('B', {});
for f = frames([5 20 77])
    frameCh2 = movCh3(:,:,f);
    [mask, blob, imgFilt] = blobMask(frameCh2, thresh3);
    statsMask = regionprops(mask, 'Area', 'PixelIdxList');
    maskMain = zeros(size(mask));
    maskMain(statsMask(find([statsMask.Area] == max([statsMask.Area]))).PixelIdxList) = 1;
    B = bwboundaries(maskMain);
    maskFrames(:,:,f) = maskMain;
    bFrames(f).B = B{1};
    
    pctl = prctile(frameCh2(maskMain == 1), 70);
    pctlMask = (frameCh2 > pctl).*maskMain;
    
    figure
    subplot(1,2,1)
    imshowpair(imadjust(movCh2(:,:,f), scaleCh2), maskMain)
    hold on
    boundary = B{1};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', .5)
    hold off
    subplot(1,2,2)
    imshow(maskMain+pctlMask,[0 2])
    
end
