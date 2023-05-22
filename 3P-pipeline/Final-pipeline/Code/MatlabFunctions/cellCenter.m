function [centrCell] = cellCenter(maskCell)
        % This function finds the centroid of the cell based on the mask of
        % the cell. If the centroid falls outside of the mask, the nearest
        % point on the mask is used instead.
        
        statsCell = regionprops(maskCell, 'Centroid', 'PixelIdxList');
        [cellRow, cellCol] = ind2sub(size(maskCell), statsCell.PixelIdxList);
        cellCoord(:,1) = cellCol; % X
        cellCoord(:,2) = cellRow; % Y
        centroidCoord = [round(statsCell.Centroid(1,1)), round(statsCell.Centroid(1,2))];
        if maskCell(round(statsCell.Centroid(1,2)), round(statsCell.Centroid(1,1))) == 0
            centrCellIdx = dsearchn(cellCoord, centroidCoord);
            centrCell = [cellCoord(centrCellIdx,1),cellCoord(centrCellIdx,2)];
        else
            centrCell = centroidCoord;
        end
end