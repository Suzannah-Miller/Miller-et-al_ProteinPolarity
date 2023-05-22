function [maskRotated, centroidsRotated, smoothProfile, intensityHalves, profileScore, brightPixHalf] = profile(cellBg, maskCell, maskPctlClose, centrCell, centrWramp, angle, option)
    % This function calculates the profile score, fold change in mean
    % intensity for halves of the cell, and the proportion of bright
    % pixels in the half of the cell.

    % Rotate cell so that axis between cell centroid and WRAMP centroid
    % lies along x-axis
    imgRotated = imrotate(cellBg, -angle);
    centroids = zeros(size(cellBg));
    centroids(round(centrCell(1,2)),round(centrCell(1,1))) = 1;
    centroids(round(centrWramp(2)),round(centrWramp(1))) = 1;
    centroids = imdilate(centroids, strel('disk', 2));
    maskRotated = imrotate(maskCell, -angle);
    centroidsRotated = imrotate(centroids, -angle);
    
    % Sum columns in image and find first and last column in rotated cell.
    startCol = find(sum(maskRotated) > 0, 1, 'first');
    endCol = find(sum(maskRotated) > 0, 1, 'last');
    halfCol = startCol + round(0.5*(endCol - startCol));
    imgRotatedNaN = imgRotated;
    imgRotatedNaN(maskRotated == 0) = NaN;
    if option == "WRAMP"
        intensityHalves = mean(imgRotatedNaN(:,halfCol:endCol),'all','omitnan')/mean(imgRotatedNaN(:,startCol:(halfCol-1)),'all','omitnan');
    elseif option == "Major"
        half1 = mean(imgRotatedNaN(:,halfCol:endCol),'all','omitnan');
        half2 = mean(imgRotatedNaN(:,startCol:(halfCol-1)),'all','omitnan');
        if half1 > half2
            intensityHalves = half1/half2;
        else
            intensityHalves = half2/half1;
        end
    end
    % Calculate intensity profiles (average values along WRAMP axis)
    profile = mean(imgRotatedNaN, 1, 'omitnan');
    smoothProfile = movmean(profile, 10);
    profileScore = (max(smoothProfile(halfCol:endCol)) - min(smoothProfile(startCol:(halfCol-1))))/(max(smoothProfile) - min(smoothProfile));
    if profileScore < 0 %>% occurs rarely for major axis due to implementation with regards to defining the major axis angle and not accounting for which side is brighter
        profileScore = 0;
    end
    % Calculate the proportion of bright pixels in halfs of the cell
    maskPctlCloseRotated = imrotate(maskPctlClose, -angle);
    brightPixHalf = sum(maskPctlCloseRotated(:,halfCol:endCol),'all')/sum(maskPctlCloseRotated(:));
end