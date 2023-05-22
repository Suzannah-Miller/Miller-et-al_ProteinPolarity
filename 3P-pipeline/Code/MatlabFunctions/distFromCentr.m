function [distTable, centrWramp] = distFromCentr(maskCombined, centrCell, option)
    % Returns table of distance measurements, including distCentr, 
    % the distance between the rounded [X,Y] 
    % coordinates for the centroid of a segmented cell (centrCell) 
    % and the rounded coordinates for the centroid of the segmented
    % brightest object (centrWramp).
    % Also returns the angle of a line
    % drawn between the center of the cell and center of the brightest object
    % (wrampAngle). Additionally, 
    % this function calculates the distance from the cell center 
    % to every point in the brightest object and returns descriptive statistics. 

    % Note: if the centroid of the cell or wramp structure falls 
    % outside of the object, the nearest point on the mask is used instead.

    if option == "WRAMP"
        maskWramp = zeros(size(maskCombined));
        maskWramp(maskCombined == 3) = 1;
    elseif option == "Second"
        maskWramp = zeros(size(maskCombined));
        maskWramp(maskCombined == 2) = 1;
    end
    clear goodY goodX x y 
    
    if sum(maskWramp, 'all') > 0
        statsWramp = regionprops(maskWramp, 'Area', 'Centroid', 'PixelIdxList');
        [wrampRow, wrampCol] = ind2sub(size(maskWramp), statsWramp.PixelIdxList);
        wrampCoord(:,1) = wrampCol; % X
        wrampCoord(:,2) = wrampRow; % Y

        % Find centrWramp (rounded centroid of maskObject or pixel in
        % maskObject closest to centroid if centroid is outside of WRAMP
        % structure mask)
        centroidCoord = [round(statsWramp.Centroid(1,1)), round(statsWramp.Centroid(1,2))];
        if maskWramp(round(statsWramp.Centroid(1,2)), round(statsWramp.Centroid(1,1))) == 0
            wrampCellIdx = dsearchn(wrampCoord, centroidCoord);
            centrWramp = [wrampCoord(wrampCellIdx,1),wrampCoord(wrampCellIdx,2)];
        else
            centrWramp = centroidCoord;
        end

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

            distCentrs = sqrt((centrWramp(1)-centrCell(1))^2 + (centrWramp(2)-centrCell(2))^2);
        end

        dist2pixel = double.empty(length(wrampRow),0);
        for i = 1:length(wrampRow)
            if wrampRow(i) == centrCell(2) && wrampCol(i) == centrCell(1)
                dist2pixel(i) = 1;
            else
                dist2pixel(i) = sqrt((wrampCol(i)-centrCell(1))^2 + (wrampRow(i)-centrCell(2))^2);
            end
        end

        distTable = table;
        distTable.minDistCentr2Wramp = min(dist2pixel);
        distTable.maxDistCentr2Wramp = max(dist2pixel);
        distTable.medianDistCentr2Wramp = median(dist2pixel);
        distTable.firstQuartileDistCentr2Wramp = prctile(dist2pixel, 25);
        distTable.thirdQuartileDistCentr2Wramp = prctile(dist2pixel, 75);
        distTable.distCentrs = distCentrs;
        distTable.wrampAngle = wrampAngle;
    else
        distTable = table('Size', [1 7], 'VariableTypes', repmat("doubleNaN", 1, 7), 'VariableNames', {'minDistCentr2Wramp','maxDistCentr2Wramp','medianDistCentr2Wramp','firstQuartileDistCentr2Wramp','thirdQuartileDistCentr2Wramp','distCentrs','wrampAngle'});
        centrWramp = NaN;
    end

    
%     figure
%     imshow(maskCombined, [0 3])
%     hold on
%     plot(centrCell(1), centrCell(2), '.','color','blue')
%     plot(centrWramp(1), centrWramp(2), '.','color','red')
%     line([centrWramp(1) centrCell(1)], [centrWramp(2) centrCell(2)], 'color', 'green')
%     hold off
end