function [mask, blob, imgFilt] = blobMask(img, thresh)
    blob = imgaussfilt(img, 4);
    blob = imdilate(blob, strel('diamond', 10));

    imgFilt = imsubtract(img, 0.15.*blob);
    imgFilt = medfilt2(imgFilt, [5 5]);

    % Threshold
    mask = imgFilt > thresh;
end