function [colocTable] = colocBrightest3ch2(maskPctlMain, maskPctlOther, maskCombinedMain, maskCombinedOther, cellBgMain, cellBgOther, threshMain, threshOther, centrCell, wrampAngleCheck)
    % One correction from 20220331: manders for the whole cell was using
    % the cell mask instead of maskPctl so all values returned were 1. This
    % was corrected in the three-channel version but not the two-channel
    % function.

    % This function calculates the colocalization features for a
    % two-channel image
    
    % Percent area overlap (relative to MCAM or F-actin) of percentile mask
    percOverlapOther2Main = sum((maskPctlMain.*maskPctlOther), 'all')/sum(maskPctlMain, 'all')*100;
    percOverlapThird2Main = NaN;
    percOverlapOther2Third = NaN;
    percOverlapMain2Other = sum((maskPctlMain.*maskPctlOther), 'all')/sum(maskPctlOther, 'all')*100;
    percOverlapMain2Third = NaN;
    percOverlapThird2Other = NaN;
    percOverlapTriple = NaN;
    
    maskCell = maskCombinedMain > 0;
    
    % Percent area overlap (relative to MCAM or F-actin) of WRAMP
    % object(brightest integrated intensity)
    maskWrampMain = zeros(size(maskCombinedMain));
    maskWrampMain(maskCombinedMain == 3) = 1;

    maskWrampOther = zeros(size(maskCombinedOther));
    maskWrampOther(maskCombinedOther == 3) = 1;
    
    percOverlapBrightestThird2Main = NaN;
    percOverlapBrightestMain2Third = NaN;
    
    percOverlapBrightestOther2Main = sum(maskWrampMain.*maskWrampOther,'all')/sum(maskWrampMain,'all')*100;
    percOverlapBrightestMain2Other = sum(maskWrampMain.*maskWrampOther,'all')/sum(maskWrampOther,'all')*100;
    
    percOverlapBrightestOther2Third = NaN;
    percOverlapBrightestThird2Other = NaN;
    
    % Percent area of percentile mask for other channel at WRAMP structure
    thirdAtWramp = NaN;
    otherAtWramp = sum((maskWrampMain.*maskPctlOther), 'all')/sum(maskPctlOther, 'all')*100;
    otherAtWrampThird = NaN;
    
    % Pearson coefficient at WRAMP structure (MCAM segmentation)
    pearsonThirdWramp = NaN;
    pearsonOtherWramp = pearson(cellBgMain(maskWrampMain == 1), cellBgOther(maskWrampMain == 1));
    pearsonOtherThirdWramp = NaN;
    
    % Manders coefficients for brightest object (should be the same as
    % percOverlapBrightest)
    mandersMainThirdBrightest1 = NaN;
    mandersMainThirdBrightest2 = NaN;
    [mandersMainOtherBrightest1, mandersMainOtherBrightest2] = manders(cellBgMain, cellBgOther, maskWrampMain, maskWrampOther, threshMain, threshOther);
    mandersThirdOtherBrightest1 = NaN;
    mandersThirdOtherBrightest2 = NaN;
    
    % Pearson/Manders coefficient for whole cell
    pearsonMainThird = NaN;
    pearsonMainOther = pearson(cellBgMain(maskCell == 1), cellBgOther(maskCell == 1));
    pearsonOtherThird = NaN;
    mandersMainThird1 = NaN;
    mandersMainThird2 = NaN;
    [mandersMainOther1, mandersMainOther2] = manders(cellBgMain, cellBgOther, maskPctlMain, maskPctlOther, threshMain, threshOther);
    mandersThirdOther1 = NaN;
    mandersThirdOther2 = NaN;
    
    % Pearson coefficient for cell excluding WRAMP structure
    pearsonMainThirdNotWramp = NaN;
    pearsonMainOtherNotWramp = pearson(cellBgMain((maskCombinedMain > 0) &(maskCombinedMain < 3)), cellBgOther((maskCombinedMain > 0) &(maskCombinedMain < 3)));
    
    % Fold-change mean intensity at WRAMP structure over rest of cell
    thirdWrampOverNot = NaN;
    otherWrampOverNot = mean(cellBgOther(maskWrampMain == 1), 'all')/mean(cellBgOther((maskCombinedMain > 0) &(maskCombinedMain < 3)), 'all');
    
    % Distance between centroids or weighted centroids (value is negative
    % if other channel's centroid is closer to the cell center than MCAM)
    distCentrThird2Main = NaN;
    distWeightedThird2Main = NaN;
    sameSideCenterThird2Main = NaN;
    sameSideWeightedThird2Main = NaN;
    [distCentrOther2Main, distWeightedOther2Main, sameSideCenterOther2Main, sameSideWeightedOther2Main] = dist2centers(maskWrampMain, maskWrampOther, maskCell, cellBgMain, cellBgOther);
    distCentrOther2Third = NaN;
    distWeightedOther2Third = NaN;
    sameSideCenterOther2Third = NaN;
    sameSideWeightedOther2Third = NaN;
 
    colocTable = table(percOverlapOther2Main, percOverlapThird2Main, percOverlapOther2Third, percOverlapMain2Other, percOverlapMain2Third, percOverlapThird2Other, percOverlapTriple, percOverlapBrightestThird2Main, percOverlapBrightestMain2Third, percOverlapBrightestOther2Main, percOverlapBrightestMain2Other, percOverlapBrightestOther2Third, percOverlapBrightestThird2Other, thirdAtWramp, otherAtWramp, otherAtWrampThird, pearsonThirdWramp, pearsonOtherWramp, pearsonOtherThirdWramp, mandersMainThirdBrightest1, mandersMainThirdBrightest2, mandersMainOtherBrightest1, mandersMainOtherBrightest2, mandersThirdOtherBrightest1, mandersThirdOtherBrightest2, pearsonMainThird, pearsonMainOther, pearsonOtherThird, mandersMainThird1, mandersMainThird2, mandersMainOther1, mandersMainOther2, mandersThirdOther1, mandersThirdOther2, pearsonMainThirdNotWramp, pearsonMainOtherNotWramp, thirdWrampOverNot, otherWrampOverNot, distCentrThird2Main, distCentrOther2Main, distCentrOther2Third, distWeightedThird2Main, distWeightedOther2Main, distWeightedOther2Third, sameSideCenterThird2Main, sameSideWeightedThird2Main, sameSideCenterOther2Main, sameSideWeightedOther2Main, sameSideCenterOther2Third, sameSideWeightedOther2Third);
    
    function [mandersCoeff1, mandersCoeff2] = manders(img1, img2, mask1, mask2, thresh1, thresh2)
        intensitiesThresh1 = (img1 - thresh1).*mask1; % Only consider pixels above threshold within object of interest
        intensitiesThresh1(intensitiesThresh1 < 0) = 0;
        intensitiesThresh2 = (img2 - thresh2).*mask2;
        intensitiesThresh2(intensitiesThresh2 < 0) = 0;
        
        mandersCoeff1 = sum(intensitiesThresh1(mask2 == 1), 'all')/sum(intensitiesThresh1(mask1 == 1), 'all');
        mandersCoeff2 = sum(intensitiesThresh2(mask1 == 1), 'all')/sum(intensitiesThresh2(mask2 == 1), 'all');
    end
    function [pearsonCoeff] = pearson(intensities1, intensities2)
        pears = corrcoef(intensities1, intensities2);
        if length(intensities1) == 1
            pearsonCoeff = 1;
        else
            pearsonCoeff = pears(1,2);
        end
    end

    function [distCenterWramp, distWeighted, sameSideCenter, sameSideWeighted] = dist2centers(maskWramp, maskWrampOther, maskCell, img1, img2)
        combMask = uint8(maskCell) + uint8(maskWramp).*2 + uint8(maskWrampOther);
        
        centrWramp = cellCenter(maskWramp);

        if centrWramp(1) == centrCell(1) && centrWramp(2) == centrCell(2)
            % Values returned if center of cell and center of wramp candidate 
            % are the same.
            wrampAngle = 0;
            distCentrs = 1;
        else % If centrWramp is not equal to centrCell (most cases)
            % Find angle of line drawn between centrCell and centrWramp     
            m = (centrWramp(2)-centrCell(2))/(centrWramp(1)-centrCell(1)); % m = slope

            if m == 0
                if centrWramp(1) > centrCell(1)
                    wrampAngle = 0;
                else
                    wrampAngle = 180;
                end
            elseif abs(m) == Inf
                if centrCell(2) > centrWramp(2)
                    wrampAngle = 90;
                else
                    wrampAngle = 270;
                end
            elseif -m > 0
                if centrCell(2) > centrWramp(2)
                    wrampAngle = atan(-m)*(180/pi);
                else
                    wrampAngle = (atan(-m)*(180/pi)) + 180;
                end
            else
                if centrWramp(1) > centrCell(1)
                    wrampAngle = (atan(-m)*(180/pi)) + 360;
                else
                    wrampAngle = (atan(-m)*(180/pi)) + 180;
                end
            end
        end
        
        if all(maskWramp == maskWrampMain)
            if wrampAngle ~= wrampAngleCheck
                error("wrampAngle doesn't match input")
            end
        end
        
        combRotated = imrotate(combMask, -wrampAngle);
        cellRotated = imdilate(imrotate(maskCell, -wrampAngle), strel('diamond',1));
%         if sum(wrampRotated, 'all') == 0
%             wrampRotated = imrotate(imdilate(maskWramp, strel('diamond',1)), -wrampAngle);
%         end
        
        maskWrampRotated = zeros(size(combRotated));
        maskWrampOtherRotated = zeros(size(combRotated));
        maskWrampRotated(combRotated == 3) = 1;
        maskWrampRotated(combRotated == 4) = 1;
        maskWrampOtherRotated(combRotated == 2) = 1;
        maskWrampOtherRotated(combRotated == 4) = 1;
       
        centrWrampRotated = cellCenter(maskWrampRotated);
        centrOtherRotated = cellCenter(maskWrampOtherRotated);
        centrCellRotated = cellCenter(cellRotated);
        
%         if centrWrampRotated(1) < (centrCellRotated(1) - 1)
%             error("Error - wrampAngle: wramp to left of center")
%         end
        distCenterWramp = sqrt((centrWrampRotated(1)-centrOtherRotated(1))^2 + (centrWrampRotated(2)-centrOtherRotated(2))^2);
        if abs(distCenterWramp) < 1
            distCenterWramp = 1;
        end
        if centrOtherRotated(1) < centrWrampRotated(1) % Centroid of other object is to the left of WRAMP
            distCenterWramp = -distCenterWramp;  
        end
        if centrOtherRotated(1) < centrCellRotated(1) && centrCellRotated(1) < centrWrampRotated(1)
            sameSideCenter = 0;
        elseif centrWrampRotated(1) <= centrCellRotated(1) && centrOtherRotated(1) < centrWrampRotated(1)
            sameSideCenter = 0;
        else
            sameSideCenter = 1;
        end   
        img1Rotated = imrotate(img1, -wrampAngle);
        img2Rotated = imrotate(img2, -wrampAngle);
        
        statsMain = regionprops(maskWrampRotated, img1Rotated, 'WeightedCentroid');
        centrMassWramp = [statsMain.WeightedCentroid];
        statsOther = regionprops(maskWrampOtherRotated, img2Rotated, 'WeightedCentroid');
        centrMassOther = [statsOther.WeightedCentroid];
        distWeighted = sqrt((centrMassWramp(1)-centrMassOther(1))^2 + (centrMassWramp(2)-centrMassOther(2))^2);
        if distWeighted < 1
            distWeighted = 1;
        end
        if centrMassOther(1) < centrMassWramp(1)
            distWeighted = -distWeighted;
        end
        if centrMassOther(1) < centrCellRotated(1) && centrCellRotated(1) < centrWrampRotated(1)
            sameSideWeighted = 0;
        elseif centrWrampRotated(1) <= centrCellRotated(1) && centrMassOther(1) < centrWrampRotated(1)
            sameSideWeighted = 0;
        else
            sameSideWeighted = 1;
        end
%         figure
%         imshow(combRotated, [0 4])
%         hold on
%         plot(centrCellRotated(1), centrCellRotated(2), '.', 'Color','blue')
%         plot(centrWrampRotated(1), centrWrampRotated(2), '.', 'Color', 'red')
%         plot(centrOtherRotated(1), centrOtherRotated(2), '.', 'Color', 'green')
%         line([centrCellRotated(1) centrWrampRotated(1)], [centrCellRotated(2) centrWrampRotated(2)], 'Color', 'black')
%         
%         figure
%         imshow(combRotated, [0 4])
%         hold on
%         plot(centrCellRotated(1), centrCellRotated(2), '.', 'Color','blue')
%         plot(centrMassWramp(1), centrMassWramp(2), '.', 'Color', 'red')
%         plot(centrMassOther(1), centrMassOther(2), '.', 'Color', 'green')
%         line([centrCellRotated(1) centrMassWramp(1)], [centrCellRotated(2) centrMassWramp(2)], 'Color', 'black')
    end
end