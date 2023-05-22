%% Introduction

% This script is used to generate masks of the images using the "blobMask"
% function and save final segmentation images. 

% The first chunk of this script reads the metadata from images, saves
% metadata to a .csv file, and checks for constant exposure, etc.

% The global threshold for the blobMask function is then set manually and 
% adjusted on a per image basis. This code shows the settings for the first
% replicate of data shown in Figures 7 and 9 as an example.

% These masks are saved as a binary image and as an overlay on the 
% fluorescence image.
% The overlay is then manually corrected and individual cells are
% segmented in ImageJ.

% The last chunk of this code is then run to use the corrected masks to
% save .tifs of the label matrices of the segmented image.

%%
clear
close all

addpath '/Applications/MATLAB_R2020a.app/toolbox/bfmatlab' % BioFormats toolbox
addpath '/.../MatlabFunctions' % Folder with custom functions

% Image directories
imgDir1 = "/.../CrisprClones/ImageData/NikonWidefield/Pilot/20220907"; % Images analyzed using machine learning.
imgDir2 = ""; % Enter additional directories if images are contained in more than one folder.
imgDir3 = ""; % Enter additional directories if images are contained in more than one folder, etc.

imgDir = [imgDir1 imgDir2 imgDir3];

% Path to directory where analyses will be saved.
directorySave = '/.../CrisprClones/Analysis/MatlabAnalysis/';
% Images to exclude from analysis
notQuantifiedImg = ["20220907_KKGKclone4_MCAM488-120_MSN594-100_Ph405-150_60X_001"]; % Longer MCAM exposure than other images

cd(directorySave)
mkdir Masks
%% Read metadata 
% Read metadata from .nd2 file and save .csv with relevant metadata for each image.
for d = 1:numel(imgDir)
    cd(imgDir(d))
    
    imgInfo = dir('*.nd2*');

    imnames = string({imgInfo.name});
    
    if isempty(notQuantifiedImg) ~= 1
        imnames(contains(imnames, notQuantifiedImg)) = [];
    end
    metaTable = table('Size',[length(imnames) 11],'VariableTypes',{'string','string','double','string','double','string','double','string','double','double','double'}, 'VariableNames', {'Image','DateAcquired','ExposureCh1','PowerCh1','ExposureCh2','PowerCh2','ExposureCh3','PowerCh3', 'DAPI', 'FITC', 'TxRedCy5'});
    
    for imnum = 1:length(imnames) % [2 12]
        clearvars imgbf imgcell img

        imgbf = bfopen(char(imnames(imnum)));
        imgcell = imgbf{1,1};
        for j = 1:size(imgcell, 1) %j = number of frames; 
            img(:,:,j) = im2double(imgcell{j,1});
        end
        metaTable.Image(imnum) = imnames(imnum);
        origMeta = imgbf{1,2};
        omeMeta = imgbf{1,4};
        
        allKeys = arrayfun(@char, origMeta.keySet.toArray, 'UniformOutput', false);
        allValues = cellfun(@(x) origMeta.get(x), allKeys, 'UniformOutput', false);
        dateKey = find(contains(string(allKeys), "TextInfoItem_9"));
        metaTable.DateAcquired(imnum) = extractAfter(string(allKeys{dateKey}), "TextInfoItem_9");
        E1 = origMeta.get('Global Exposure time (text) #1');
        E2 = str2double(extractBetween(string(omeMeta.getPlaneExposureTime(0,0)), "value[", "], unit[s]"));
        if E1 == E2
            metaTable.ExposureCh1(imnum) = E1;
        else
            error("Channel not identified correctly in metadata")
        end
        metaTable.PowerCh1(imnum) = string(origMeta.get('Global X-Cite XYLIS, Illuminator(XYLIS) Iris intensity #1'));
        E1 = origMeta.get('Global Exposure time (text) #2');
        E2 = str2double(extractBetween(string(omeMeta.getPlaneExposureTime(0,1)), "value[", "], unit[s]"));
        if E1 == E2
            metaTable.ExposureCh2(imnum) = E1;
        else
            error("Channel not identified correctly in metadata")
        end
        metaTable.PowerCh2(imnum) = string(origMeta.get('Global X-Cite XYLIS, Illuminator(XYLIS) Iris intensity #2'));
        if j == 3
            E1 = origMeta.get('Global Exposure time (text) #3');
            E2 = str2double(extractBetween(string(omeMeta.getPlaneExposureTime(0,2)), "value[", "], unit[s]"));
            if E1 == E2
                metaTable.ExposureCh3(imnum) = E1;
            else
                error("Channel not identified correctly in metadata")
            end
            metaTable.PowerCh3(imnum) = string(origMeta.get('Global X-Cite XYLIS, Illuminator(XYLIS) Iris intensity #3'));
        end
        channelNames = [string(origMeta.get('Global Nikon Ti, FilterChanger(Turret1) #1')), string(origMeta.get('Global Nikon Ti, FilterChanger(Turret1) #2')), string(origMeta.get('Global Nikon Ti, FilterChanger(Turret1) #3'))];
        if j == 3
            metaTable.DAPI(imnum) = find(contains(channelNames, "DAPI"));
        else 
            metaTable.DAPI(imnum) = 0;
        end
        metaTable.FITC(imnum) = find(contains(channelNames, "FITC"));
        metaTable.TxRedCy5(imnum) = find(contains(channelNames, ["Texas", "Cy5"]));
    end
    if d == 1
        metaTableAll = metaTable;
    else
        metaTableAll = vertcat(metaTableAll, metaTable);
    end
end
writetable(metaTableAll, [directorySave, 'BeginAnnotationFilePilot.csv'])

for c = 3:8
    if height((unique(metaTableAll(:,c)))) > 1
        disp(['Warning: Exposure not constant for all images. Column ', num2str(c)])
    end
end
for c = 9:11
    if height((unique(metaTableAll(:,c)))) > 1
        disp(['Warning: Channel order is not the same for all images. Column ', num2str(c)])
    end
end

for r = 1:height(metaTableAll)
    fitcExp(r) = metaTableAll{r,['ExposureCh', num2str(metaTableAll.FITC(r))]};
    txRedCy5Exp(r) = metaTableAll{r,['ExposureCh', num2str(metaTableAll.TxRedCy5(r))]};
end

tempTable = metaTableAll;
noDapiIdx = tempTable.DAPI == 0;
tempTable(noDapiIdx,:) = [];
for r = 1:height(tempTable)
    dapiExp(r) = tempTable{r,['ExposureCh', num2str(tempTable.DAPI(r))]};
end

if length(unique(dapiExp)) > 1
    disp('Warning: DAPI exposure is not constant for all images')
end

if length(unique(fitcExp)) > 1
    disp('Warning: FITC exposure is not constant for all images')
end

if length(unique(txRedCy5Exp)) > 1
    disp('Warning: Texas Red exposure is not constant for all images')
end

if length(unique(metaTableAll.Image)) ~= height(metaTableAll)
    disp('Warning: Image names are not unique for this experiment')
end

%% Generate masks for manual correction
% Combines masks from two image channels

% Only ran d = 1 for now
annotation = readtable([directorySave, 'BeginAnnotationFilePilot.csv']);
for d = 1%:numel(imgDir)
    cd(imgDir(d))
    
    imgInfo = dir('*.nd2*');

    imnames = string({imgInfo.name});
    
    imgInfo = dir('*.nd2*');

    imnames = string({imgInfo.name});
    if isempty(notQuantifiedImg) ~= 1
        imnames(contains(imnames, notQuantifiedImg)) = [];
    end
    
    for imnum = 1:length(imnames)
        clearvars imgbf imgcell img imgMcam imgFactin imgMcamAdj imgFactinAdj
        annotationTemp = annotation(contains(annotation.Image, imnames(imnum)),:);
        
        imgbf = bfopen(char(imnames(imnum)));
        imgcell = imgbf{1,1};
        for j = 1:size(imgcell, 1) %j = number of frames; 
            img(:,:,j) = im2double(imgcell{j,1});
        end
        
        % Use this code to set the threshold for blobMask and adjust as
        % needed for different images.
        if contains(annotationTemp.Image, "MCAM488")
            mcamChannel = annotationTemp.FITC;
            threshM = 0.003;
            threshF = 0.006;
            if contains(annotationTemp.Image, ["20220907_EricR23_MCAM488-100_MSN594-100_Ph405-150_60X_S1_003", "20220907_EricR23_MCAM488-100_MSN594-100_Ph405-150_60X_S1_004", "20220907_EricR23_MCAM488-100_MSN594-100_Ph405-150_60X_S1_005", "20220907_EricR23_MCAM488-100_MSN594-100_Ph405-150_60X_S1_006"])
                threshM = 0.004;
                threshF = 0.009;
            end
            if contains(annotationTemp.Image, ["20220907_KKGKclone3_MCAM488-100_MSN594-100_Ph405-150_60X_S1"])
                threshM = 0.00325;
                threshF = 0.008;
            end
            if contains(annotationTemp.Image, ["20220907_KKGKclone3_MCAM488-100_MSN594-100_Ph405-150_60X_S2"])
                threshM = 0.00375;
                threshF = 0.00875;
            end
            if contains(annotationTemp.Image, ["20220907_KKGKclone4_MCAM488-100_MSN594-100_Ph405-150_60X_S1"])
                threshM = 0.00325;
                threshF = 0.007;
            end
            if contains(annotationTemp.Image, ["20220907_KKGKclone4_MCAM488-100_MSN594-100_Ph405-150_60X_S2"])
                if contains(annotationTemp.Image, ["20220907_KKGKclone4_MCAM488-100_MSN594-100_Ph405-150_60X_S2_001_5"])
                    threshM = 0.003;
                    threshF = 0.00575;
                elseif contains(annotationTemp.Image, ["20220907_KKGKclone4_MCAM488-100_MSN594-100_Ph405-150_60X_S2_001_1", "20220907_KKGKclone4_MCAM488-100_MSN594-100_Ph405-150_60X_S2_001_2"]) 
                    threshM = 0.002;
                    threshF = 0.00425;
                else
                    threshM = 0.0022;
                    threshF = 0.004;
                end
            end
        elseif contains(annotationTemp.Image, "EZRM488")
            mcamChannel = annotationTemp.FITC;
            threshM = 0.00375;
            threshF = 0.006;
        elseif contains(annotationTemp.Image, "MSNM488")
            mcamChannel = annotationTemp.TxRedCy5; % Use MCAM (rabbit) channel
            threshM = 0.00175;
            threshF = 0.006;
        elseif contains(annotationTemp.Image, "MSN594") % Control (MSN only)
            mcamChannel = annotationTemp.TxRedCy5;
            threshM = 0.003;
            threshF = 0.006;
        else
            mcamChannel = annotationTemp.DAPI; % Control (Phalloidin only)
            threshM = 0.006;
            threshF = 0.006;
        end
            
        if contains(annotationTemp.Image, "Ph405")
            if contains(annotationTemp.Image, "EricR23_488-100ms_594-100ms_Ph405")
                secondMask = 0;
            else
                secondMask = 1;
            end
        else
            secondMask = 0;
        end
        factinChannel = annotationTemp.DAPI;
        
        imgMcam = img(:,:,mcamChannel);
        imgMcamAdj = imadjust(imgMcam, [0 0.1]);
       
        % Generate mask for MCAM channel using blobMask
        [maskM, blobM, imgFiltM] = blobMask(imgMcam, threshM);
        
%         figure % Uncomment this code to test different thresholds
%         imshow([imgMcamAdj imadjust(blobM, [0 0.04]) imadjust(imgFiltM, [0 0.04])])
%         testThresh = linspace(0.002,0.003, 5);
%         close all
%         for i = 1:length(testThresh)
%             maskMtest = imgFiltM > testThresh(i);
%             figure
%             imshowpair(imgMcamAdj, 0.3.*maskMtest, 'Scaling', 'none')
%             title(['Threshold = ', num2str(testThresh(i))])
%         end
        
        if secondMask == 1 % If using a second image channel (F-actin) to generate mask
            imgFactin = img(:,:,factinChannel);
            imgFactinAdj = imadjust(imgFactin, [0 0.05]);

%     %         figure
%     %         imshow([imgMcamAdj imgFactinAdj])
%             
            [maskF, blobF, imgFiltF] = blobMask(imgFactin, threshF);
%  
%             figure % Uncomment this code to test different thresholds
%             imshow([imgFactinAdj imadjust(blobF, [0 0.04]) imadjust(imgFiltF, [0 0.04])])
%             testThresh = linspace(0.004,.006, 5);
%             close all
%             for i = 1:length(testThresh)
%                 maskFtest = imgFiltF > testThresh(i);
%                 figure
%                 imshowpair(imgFactinAdj, 0.3.*maskFtest, 'Scaling', 'none')
%                 title(['Threshold = ', num2str(testThresh(i))])
%             end
        else
            maskF = zeros(size(maskM));
        end
        maskCombined = (maskM + maskF) > 0;
%         figure
%         imshowpair(maskM, maskF)
%         figure
%         imshowpair(imgMcamAdj, 0.3.*maskCombined, 'Scaling', 'none')

        imgOverlay = (imadjust(imgMcam, [0 prctile(imgMcam(maskCombined), 75)]) + 0.5.*maskCombined)./2;
        % Save mask and overlay image.
        imwrite(maskCombined, [directorySave, '/Masks/', char(extractBefore(imnames(imnum),'.nd2')), '.tif'], 'Compression', 'none');
        imwrite(uint8(imgOverlay.*255), [directorySave, '/Masks/', char(extractBefore(imnames(imnum),'.nd2')), '_overlay.tif'], 'Compression', 'none');

    end
end

%% Segmentation
cd(directorySave)

mkdir ManualSegmentation
mkdir SegmentationCorrectionFigures

cd([directorySave, 'CorrectedMasks'])
correctedOverlayInfo = dir('*overlay*');
overlaynames = string({correctedOverlayInfo.name});

%imgNoCells = [];

cd([directorySave, 'Masks'])
maskInfo = dir('*.tif*');
masknames = string({maskInfo.name});
%masknames(contains(masknames, imgNoCells)) = [];

for imnum = 1:length(correctedOverlayInfo)
    cd([directorySave, 'Masks'])
    maskIdx = find(matches(masknames, [char(extractBefore(overlaynames(imnum), "_overlay")), '.tif']));

    originalMask = imread(masknames(maskIdx(1)));

    cd([directorySave, 'CorrectedMasks'])
    overlayCorrected = im2double(imread(overlaynames(imnum)));
%     overlayCorrected(1:8070,8054) = 1;
%     overlayCorrected(8070,1:8054) = 1;

    rawMask = originalMask;
    rawMask(overlayCorrected == 0) = 0;
    rawMask(overlayCorrected == 1) = 1;

    rawMask = bwareaopen(rawMask, 11257); %40X = 0.1625 um/pixel. 60X = 0.1083 um/pixel. 5000 pixels (40X) = 11685 pixels (60X)     

    maskComp = bwareaopen(imcomplement(rawMask), 8);
    statsComp = regionprops(maskComp, 'PixelIdxList', 'Eccentricity', 'Orientation', 'MinFeretProperties');
    angleIdx = find(abs([statsComp.Orientation]) > 89);
    lineIdx = find([statsComp.Eccentricity] > 0.97);
    widthIdx = find([statsComp.MinFeretDiameter] < 10);
    badCellIdx = intersect(intersect(angleIdx, lineIdx),widthIdx);
    lineMask = zeros(size(rawMask));
    for i = 1:length(badCellIdx)
        lineMask(statsComp(badCellIdx(i)).PixelIdxList)= 1;
    end
    lineMask = imbinarize(lineMask);

    maskConn = bwconncomp(imfill(rawMask, 'holes'));
    maskCorrected = imfill(rawMask, 'holes');
    for i = 1:maskConn.NumObjects
        maskTemp = zeros(size(rawMask));
        maskTemp(maskConn.PixelIdxList{i}) = 1;
        if sum((maskTemp.*lineMask) > 0, 'all') > 0
            maskCorrected(maskConn.PixelIdxList{i}) = 0;
        end
    end
    statsMask = regionprops(maskCorrected, 'Area', 'Perimeter', 'PixelIdxList');
    tooBigIdx = find([statsMask.Area] > 2.5e5);
    maskCorrected2 = maskCorrected;
    for i = 1:length(tooBigIdx)
        maskCorrected2(statsMask(tooBigIdx(i)).PixelIdxList) = 0;    
    end

    maskCorrected2(:,8054) = 1;
    maskCorrected2(8070,:) = 1;
    correctedConn = bwconncomp(imclearborder(maskCorrected2));
    seg = uint16(zeros(size(rawMask)));
    for i = 1:correctedConn.NumObjects
        seg(correctedConn.PixelIdxList{i}) = i;
    end

    figure
    subplot(1,2,1)
    imshow(label2rgb(seg, 'jet', 'black', 'shuffle'))
    title(masknames(maskIdx(1)), 'Interpreter', 'none', 'FontSize', 6)
    subplot(1,2,2)
    imshowpair(originalMask, imclearborder(maskCorrected2))

    imwrite(im2uint16(seg), [directorySave, 'ManualSegmentation/Segmented_', char(masknames(maskIdx(1)))], 'Compression', 'none')
    f = gcf;
    exportgraphics(f, [directorySave, 'SegmentationCorrectionFigures/', extractBefore(char(masknames(maskIdx(1))), '.tif'), '.png'], 'Resolution', 1200)
%     if imnum == 1
%         nCells =  i; 
%     else
%         nCells = nCells + i;
%     end
    close
    clearvars -except imnum maskInfo masknames directorySave overlaynames correctedOverlayInfo nCells
end

% First view cleaned up masks and make additional corrections to images in CorrectedMasks folder (overwrite
% previous corrected image). Then rerun and save final segmented masks.

