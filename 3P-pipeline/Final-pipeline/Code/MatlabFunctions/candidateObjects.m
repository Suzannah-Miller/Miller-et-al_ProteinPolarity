function [numCandidates, largeObj, brightMeanObj, brightSumObj, medianObjectMeanIntensity, medianObjectSumIntensity, medianObjectArea, candidateStats] = candidateObjects(maskPctlClose, cellBg)
    % This function finds the individual objects in an image after the
    % percentile threshold has been applied and calculates features for
    % each object, using the background subtracted image for intensity
    % features.
    
    candidateStats = regionprops(maskPctlClose, 'Area', 'PixelIdxList');
    numCandidates = length([candidateStats.Area]);
    for i = 1:numCandidates
        sumIntensities(i) = sum(cellBg(candidateStats(i).PixelIdxList));
        meanIntensities(i) = mean(cellBg(candidateStats(i).PixelIdxList));
    end

    sumIntensCell = num2cell(sumIntensities);
    meanIntensCell = num2cell(meanIntensities);
    [candidateStats.sumIntensities] = sumIntensCell{:};
    [candidateStats.meanIntensities] = meanIntensCell{:};
    maxArea = max([candidateStats.Area]);
    maxMeanIntensity = max(meanIntensities);
    maxSumIntensity = max(sumIntensities);
    
    medianObjectMeanIntensity = median(meanIntensities);
    medianObjectSumIntensity = median(sumIntensities);
    medianObjectArea = median([candidateStats.Area]);

    largeObj = 0;
    brightMeanObj = 0;
    brightSumObj = 0;

    for j = 1:numCandidates
        if candidateStats(j).Area > 0.25*maxArea
            largeObj = largeObj + 1;
        end
        if meanIntensities(j) > 0.5*maxMeanIntensity
            brightMeanObj = brightMeanObj + 1;
        end
        if sumIntensities(j) > 0.25*maxSumIntensity
            brightSumObj = brightSumObj + 1;
        end
    end
end
