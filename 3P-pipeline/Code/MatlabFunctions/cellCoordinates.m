function [topRow, bottomRow, leftCol, rightCol, coord] = cellCoordinates(img, PixelIdxList)
    % Extracts rectangular image of individual cell from larger image 
    % with 1 pixel padding when possible
    [rows, columns] = ind2sub(size(img),PixelIdxList);
    coord(:,1) = columns; % X
    coord(:,2) = rows; % Y
    if min(rows) == 1
        topRow = min(rows);
    else 
        topRow = min(rows) - 1;
    end
    if max(rows) == size(img,1)
        bottomRow = max(rows);
    else
        bottomRow = max(rows) + 1;
    end
    if min(columns) == 1
        leftCol = min(columns);
    else
        leftCol = min(columns) - 1;
    end
    if max(columns) == size(img,2)
        rightCol = max(columns);
    else
        rightCol = max(columns) + 1;
    end
end