function[maskPctl, maskCombined, pctl] = cellPrctileNoClose(imgCellBg, maskCell, percentile)
    % This function applies a threshold to the background 
    % subtracted image of the cell based on the provided percentile
    % of pixel intensities to generate a mask of bright pixels in 
    % the (maskPctl). It also returns a "combined mask" of the Cell 
    % where the object with the highest integrated intensity 
    % (brightest object) has a value of 3, the second brightest object 
    % has a value of 2, and the rest of the cell has a value of one, 
    % and returns the intensity value used as the threshold (pctl).
    
    % Segment bright objects using intensity cutoff
    pctl = prctile(imgCellBg(maskCell == 1), percentile);
    maskPctl = imbinarize((imgCellBg > pctl).*maskCell);
    maskPctlConn = bwconncomp(maskPctl);
    maskPctlLabel = zeros(size(maskCell));
    for i = 1:maskPctlConn.NumObjects
        maskPctlLabel(maskPctlConn.PixelIdxList{i}) = i;
        objAreas(i) = length(maskPctlConn.PixelIdxList{i});
        sumIntensities(i) = sum(imgCellBg(maskPctlConn.PixelIdxList{i}));
    end

    candidates = 1:length(objAreas);

    % Select object with greatest integrated intensity
    brightestObjIdx = find(sumIntensities(candidates) == max(sumIntensities(candidates)));
    maskWramp = zeros(size(maskCell));
    maskWramp(maskPctlConn.PixelIdxList{candidates(brightestObjIdx)}) = 1;
    if length(candidates) > 1
        secondBrightestObjIdx = find(sumIntensities(candidates) == max(sumIntensities(candidates(candidates ~= candidates(brightestObjIdx)))));
        maskSecond = zeros(size(maskCell));
        maskSecond(maskPctlConn.PixelIdxList{candidates(secondBrightestObjIdx(1))}) = 1;
    end

    % Combine mask of whole cell and WRAMP structure; pixels in WRAMP
    % structure have a value of 3 while pixels in the second brightest object = 2, and pixels 
    % not in the WRAMP structure have a value of 1.
    maskCombined = uint8(maskWramp).*2 + uint8(maskCell);
    if length(candidates) > 1
        maskCombined = uint8(maskWramp).*2 + uint8(maskCell) + uint8(maskSecond);
    end
end


