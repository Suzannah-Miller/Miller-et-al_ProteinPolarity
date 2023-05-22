function [normProfile] = normProfileValues(maskCell, cellBg, angle)
    % Calculates the smoothed profile as it is calculated in the "profile"
    % function, and scale the values between the minimum and maximum value
    % for purposes of generating the profiles in the figure. (Not used in
    % ML classification).

    % Rotate cell so that axis between cell centroid and WRAMP centroid
    % lies along x-axis
    imgRotated = imrotate(cellBg, -angle);
    maskRotated = imrotate(maskCell, -angle);
    
    imgRotatedNaN = imgRotated;
    imgRotatedNaN(maskRotated == 0) = NaN;
    
    % Calculate intensity profiles (average values along WRAMP axis)
    profile = mean(imgRotatedNaN, 1, 'omitnan');
    smoothProfile = movmean(profile, 10);
    
    maxVal = max(smoothProfile);
    minVal = min(smoothProfile);
    normProfile = (smoothProfile - minVal)./(maxVal - minVal);
end