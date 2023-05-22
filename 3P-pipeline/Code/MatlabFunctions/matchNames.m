function[] = matchNames(imnames, masknames, segnames)
    % This function checks whether the order of images, mask images, 
    % and segmentation images is the same.
    for imnum = 1:length(imnames)
        imnameTemp = extractBefore(imnames(imnum), '.nd2');
        segnameTemp = segnames(imnum);
        masknameTemp = masknames(imnum);
        if matches(masknameTemp, [char(imnameTemp), '.tif']) == 0
            error(['Mask does not match image name (imnum = ', num2str(imnum), ').'])
        end
        if matches(segnameTemp, ['Segmented_', char(imnameTemp), '.tif']) == 0
            error(['Segmentation does not match image name (imnum = ', num2str(imnum), ').'])
        end
    end
    display('masknames match imnammes')
    display('segnames match imnames')
end