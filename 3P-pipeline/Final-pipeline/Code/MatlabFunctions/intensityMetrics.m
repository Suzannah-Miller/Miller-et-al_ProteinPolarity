function [intensityTable] = intensityMetrics(cellBg, maskCell, maskCombined, maskPctlClose)
    % Calculates mean, median, and integrated intensity and area 
    % for the cell, the brightest, and second brightest objects 
    % in the cell.
    
    % Note: "notWramp" indicates the region of the cell excluding the 
    % brightest object and "notBright" indicates the region of the cell
    % excluding the first
    % and second brightest objects. If there is not a second brightest
    % object (because all pixels above the percentile threshold were
    % contained in a single connected component), this measurement is the
    % same as "notWramp."
    
    intensityTable = table;
    intensityTable.meanIntensityCell = mean(cellBg(maskCell == 1));
    intensityTable.meanIntensityWramp = mean(cellBg(maskCombined == 3));
    intensityTable.medianIntensityCell = median(cellBg(maskCell == 1));
    intensityTable.medianIntensityWramp = median(cellBg(maskCombined == 3));
    intensityTable.sumIntensityCell = sum(cellBg(maskCell == 1));
    intensityTable.sumIntensityWramp = sum(cellBg(maskCombined == 3));
    intensityTable.areaWramp = sum(maskCombined == 3, 'all');
    intensityTable.areaBrightObj = sum(maskPctlClose(:));
    intensityTable.meanIntensityBrightObj = mean(cellBg(maskPctlClose == 1));
    intensityTable.sumIntensityBrightObj = sum(cellBg(maskPctlClose == 1));
    if sum(maskCombined == 2) == 0
        intensityTable.sumNotWramp = sum(cellBg(maskCombined == 1));
        intensityTable.meanNotWramp = mean(cellBg(maskCombined == 1));
        intensityTable.medianNotWramp = median(cellBg(maskCombined == 1));
        intensityTable.sumNotBright = intensityTable.sumNotWramp; 
        intensityTable.meanNotBright = intensityTable.meanNotWramp;
        intensityTable.medianNotBright = intensityTable.medianNotWramp;
        intensityTable.meanIntensitySecond = 0;
        intensityTable.sumIntensitySecond = 0;
        intensityTable.medianIntensitySecond = 0;
        intensityTable.areaSecond = 0;
    else
        intensityTable.sumNotWramp = sum(cellBg(maskCombined == 1 | maskCombined == 2));
        intensityTable.meanNotWramp = mean(cellBg(maskCombined == 1 | maskCombined == 2));
        intensityTable.medianNotWramp = median(cellBg(maskCombined == 1 | maskCombined == 2));
        intensityTable.sumNotBright = sum(cellBg(maskCombined == 1));
        intensityTable.meanNotBright = mean(cellBg(maskCombined == 1));
        intensityTable.medianNotBright = median(cellBg(maskCombined == 1));
        intensityTable.meanIntensitySecond = mean(cellBg(maskCombined == 2));
        intensityTable.sumIntensitySecond = sum(cellBg(maskCombined == 2));
        intensityTable.medianIntensitySecond = median(cellBg(maskCombined == 2));
        intensityTable.areaSecond = sum(maskCombined == 2, 'all');
    end
end