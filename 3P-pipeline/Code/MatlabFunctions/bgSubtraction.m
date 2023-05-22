function[imgBg, maskBg, medBg] = bgSubtraction(img, masknames, imnum)
    % Background subtraction: dilates mask of cells and takes the median
    % intensity of the background. This value is subtracted from the image.
    rawMask = imread(masknames(imnum));
    maskBg = imdilate(rawMask, strel('disk', 10));
    medBg = median(img(maskBg == 0));
    imgBg = imsubtract(img, medBg);
    imgBg(imgBg < 0) = 0;
end