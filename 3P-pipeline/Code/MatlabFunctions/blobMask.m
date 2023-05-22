function [mask, blob, imgFilt] = blobMask(img, thresh)
    % Smooth and dilate image
    blob = imgaussfilt(img, 4);
    blob = imdilate(blob, strel('diamond', 20));
    
    % Subtract from original image
    imgFilt = imsubtract(img, 0.4.*blob);
    imgFilt = medfilt2(imgFilt, [5 5]);

    % Threshold
    mask = imgFilt > thresh;
end