function [wrampEdgeDistStats] = wrampDistances(maskCombined, wrampAngle, option)
    % This function creates a mask of the perimeter of the brightest  
    % object, and measures distances along a line drawn from 
    % the cell center through each point on border of the
    % brightest object contained within +/- 45 degrees along the axis formed by the 
    % center of the cell and the center of the brightest object,
    % with the cell center as the vertex.
    % Returns descriptive statistics on the distance along this 
    % line from the cell center to maximum edge of the cell (distCentrEdge) 
    % and the outer edge of the brightest object
    % to the maximum edge of the cell (distWrampEdge).

    % This function aims to characterize the distance from the WRAMP
    % structure to the cell edge relative to the cell
    % center.
    
    if option == "WRAMP"
        maskWramp = zeros(size(maskCombined));
        maskWramp(maskCombined == 3) = 1;
    elseif option == "Second"
        maskWramp = zeros(size(maskCombined));
        maskWramp(maskCombined == 2) = 1;
    end
    
    if sum(maskWramp, 'all') > 0

        maskCell = maskCombined > 0;

        cellRotated = imdilate(imrotate(maskCell, -wrampAngle), strel('diamond',1));
        wrampRotated = imdilate(imrotate(maskWramp, -wrampAngle), strel('diamond',1));
        if sum(wrampRotated, 'all') == 0
            wrampRotated = imrotate(imdilate(maskWramp, strel('diamond',1)), -wrampAngle);
        end
 
        statsCell = regionprops(cellRotated, 'Centroid', 'PixelIdxList', 'Area');
        statsWramp = regionprops(wrampRotated, 'Centroid', 'PixelIdxList', 'Area');

        [cellRow, cellCol] = ind2sub(size(cellRotated), statsCell.PixelIdxList);
        cellCoord(:,1) = cellCol; % X
        cellCoord(:,2) = cellRow; % Y

        [wrampRow, wrampCol] = ind2sub(size(wrampRotated), statsWramp.PixelIdxList);
        wrampCoord(:,1) = wrampCol; % X
        wrampCoord(:,2) = wrampRow; % Y

        centroidCell = [round(statsCell.Centroid(1,1)), round(statsCell.Centroid(1,2))];
        if cellRotated(round(statsCell.Centroid(1,2)), round(statsCell.Centroid(1,1))) == 0
            wrampCellIdx = dsearchn(cellCoord, centroidCell);
            centrCell = [cellCoord(wrampCellIdx,1), cellCoord(wrampCellIdx,2)];
        else
            centrCell = centroidCell;
        end

        centroidWramp = [round(statsWramp.Centroid(1,1)), round(statsWramp.Centroid(1,2))];
        if wrampRotated(round(statsWramp.Centroid(1,2)), round(statsWramp.Centroid(1,1))) == 0
            wrampCellIdx = dsearchn(wrampCoord, centroidWramp);
            centrWramp = [wrampCoord(wrampCellIdx,1),wrampCoord(wrampCellIdx,2)];
        else
            centrWramp = centroidWramp;
        end

        x1 = centrCell(1,1):size(wrampRotated, 2);
        b1 = centrCell(1,2) + centrCell(1,1); % m = -1 for line at 45 degrees
        b2 = centrCell(1,2) - centrCell(1,1); % m = +1 for line at -45 degrees
        y1 = -x1 + b1; % Calculate points in line
        y1(y1 < 1) = [];
        y2 = (x1 + b2);
        y2(y2 > size(wrampRotated, 1)) = [];

        % Make inverted mask of lines at +/- 45 degrees
        blackMask = ones(size(wrampRotated));
        for l = 1:length(y1)
            blackMask(y1(l),x1(l)) = 0;
        end
        for l = 1:length(y2)
            blackMask(y2(l),x1(l)) = 0;
        end

        blackMask = imbinarize(imerode(blackMask, strel('diamond', 1))); % Thickens line
        bmStats = regionprops(blackMask, 'Area', 'PixelIdxList');
        bigIdx = find(max([bmStats.Area])); 
        blackMask(bmStats(bigIdx).PixelIdxList) = 0; % Remove area outside of [-45,45 degrees]

        borderWramp = edge(imfill(wrampRotated, 'holes')); % Find borders of brightest object contained within +/- 45 degrees
        borderWramp = borderWramp .* blackMask;
        borderCell = edge(cellRotated);

        maskBorderWramp = ((borderWramp + wrampRotated).*blackMask) > 0;
        maskBorderCell = (borderCell + cellRotated) > 0;

        % Get X,Y coordinates for pixels in borderWramp (only pixels along the
        % border of WRAMP structure contained in [-45,45 degrees])
        statsBorder = regionprops(borderWramp, 'PixelIdxList');
        pixelIdxList = statsBorder(1).PixelIdxList;
        if length(statsBorder) > 1 % True if border is not contiguous
            for i = 2:length(statsBorder)
                pixelIdxList = [pixelIdxList; statsBorder(i).PixelIdxList];
            end
        end

        [borderRow, borderCol] = ind2sub(size(maskBorderWramp), pixelIdxList);

        for j = 1:length(borderRow) % For every pixel in boarder of brightest object
            clear goodY goodX x y inCell inWramp edgeIdx wrampInEdgeIdx wrampOutEdgeIdx wrampInEdge wrampOutEdge edge

            % Find slope between center of cell and pixel j in wrampBorder
            m = (borderRow(j)-centrCell(2))/(borderCol(j)-centrCell(1)); % slope
            b = borderRow(j) - m*borderCol(j);

            % Find coordinates for pixels in the line starting from the cell
            % center and going through pixel j in wrampBorder out to edge of image.
            if borderCol(j) == centrCell(1) % Vertical line
                x = repmat(centrCell(1), 1, (length(cellRotated(:,1))-centrCell(2))*2);
                if centrCell(2) > borderRow(j)
                    y = linspace(1, centrCell(2), length(x));
                else
                    y = linspace(centrCell(2), size(cellRotated, 1), length(x));
                end
            elseif borderRow(j) == centrCell(2) % Horizontal line
                if borderCol(j) > centrCell(1)
                    x = linspace(centrCell(1), size(cellRotated,2), (size(cellRotated,2) - centrCell(1))*2);
                else
                    x = linspace(1, centrCell(1), (centrCell(1)-1)*2);
                end
                y = repmat(centrCell(2), 1, length(x));
            elseif borderCol(j) > centrCell(1) % WRAMP structure to the right of cell center
                if abs(m) > 1
                    x = linspace(centrCell(1),length(cellRotated(1,:)), (length(cellRotated(1,:))-centrCell(1))*abs(m)*2);
                else
                    x = linspace(centrCell(1),length(cellRotated(1,:)), (length(cellRotated(1,:))-centrCell(1))*2);
                end
                y = m.*x + b; 
            else % WRAMP structure to the left of cell center
                error('Point to left of cell center');
            end
            x = round(x);
            y = round(y); 

            goodYidx = find((y > 0) & (y <= size(cellRotated, 1))); % Y values (from equation of line) contained within image.
            if size(goodYidx,2) > 1
                goodY = y(goodYidx);
                goodX = x(goodYidx);
            else
                error('goodY not > 1 value')
            end

            for i = 1:length(goodX)
                inCell(i) = maskBorderCell(goodY(i), goodX(i));  % Coordinates of line that are within the mask will have a corresponding value of 1.
                inWramp(i) = maskBorderWramp(goodY(i), goodX(i)); 
            end

            edgeIdx = find(inCell == 1, 1, 'last'); % This is the "max" distance to the edge as it will find the farthes point along the line even if there is background in between (e.g., cell wraps around, has skinny protrusion, etc.) 
            wrampInEdgeIdx = find(inWramp == 1, 1, 'first');
            wrampOutEdgeIdx = find(inWramp == 1, 1, 'last');

            cellEdge = [goodX(edgeIdx), goodY(edgeIdx)];
            wrampOutEdge = [goodX(wrampOutEdgeIdx), goodY(wrampOutEdgeIdx)];

            % Distance from cell center to edge of cell along line through
            % pixel j
            distCentrEdge(j) = sqrt((cellEdge(1)-centrCell(1))^2 + (cellEdge(2)-centrCell(2))^2)-1;% Subtract 1 because of dilation to preserve protrusions
            if distCentrEdge(j) <= 1
                distCentrEdge(j) = 1;
            end
            % Distance from outer boarder of brightest object to edge of
            % cell along line through pixel j
            distWrampEdge(j) = sqrt((cellEdge(1)-wrampOutEdge(1))^2 + (cellEdge(2)-wrampOutEdge(2))^2); % used to use centrWramp
            if distWrampEdge(j) <= 1
                distWrampEdge(j) = 1;
            end
        end
        wrampEdgeDistStats = table;
        wrampEdgeDistStats.minDistCentrEdge = min(distCentrEdge);
        wrampEdgeDistStats.minDistWrampEdge = min(distWrampEdge);

        wrampEdgeDistStats.maxDistCentrEdge = max(distCentrEdge);
        wrampEdgeDistStats.maxDistWrampEdge = max(distWrampEdge);

        wrampEdgeDistStats.medianDistCentrEdge = median(distCentrEdge);
        wrampEdgeDistStats.medianDistWrampEdge = median(distWrampEdge);

        wrampEdgeDistStats.firstQuartileDistCentrEdge = prctile(distCentrEdge, 25);
        wrampEdgeDistStats.firstQuartileDistWrampEdge = prctile(distWrampEdge, 25);

        wrampEdgeDistStats.thirdQuartileDistCentrEdge = prctile(distCentrEdge, 75);
        wrampEdgeDistStats.thirdQuartileDistWrampEdge = prctile(distWrampEdge, 75);
    else
        wrampEdgeDistStats = table('Size', [1 10], 'VariableTypes', repmat("doubleNaN", 1, 10), 'VariableNames', {'minDistCentrEdge', 'minDistWrampEdge', 'maxDistCentrEdge', 'maxDistWrampEdge', 'medianDistCentrEdge', 'medianDistWrampEdge', 'firstQuartileDistCentrEdge', 'firstQuartileDistWrampEdge', 'thirdQuartileDistCentrEdge', 'thirdQuartileDistWrampEdge'});
    end
        
    
    
%     figure
%     fig1 = imshow(imdilate(borderCell, strel('disk',3)))
%     hold on
%     fig2 = imshow(imdilate(borderWramp, strel('disk',3)))
%     plot(centrCell(1),centrCell(2),'.', 'color', 'blue')
%     plot(goodX,goodY, 'color', 'red')
%     plot(x1(1:length(y1)), y1, 'green')
%     plot(x1(1:length(y2)), y2, 'green')
%     plot(cellEdge(1),cellEdge(2), '.','color', 'blue')
%     plot(wrampOutEdge(1), wrampOutEdge(2), '.','color','yellow')
%     plot(borderCol(j), borderRow(j), '.', 'color', 'cyan')
%     hold off
%     alpha(fig2, 0.4)
end