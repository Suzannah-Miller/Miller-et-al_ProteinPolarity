function[f1, f2] = parametersFigure(centrCell, onlyOneChannel, thirdChannel, cellMainAdjust, cellBgMain, maskCombinedMain, maskPctlMain, pctlMain, wrampProfileMain, candidateStatsMain, figNameMain, cellOtherAdjust, cellBgOther, maskCombinedOther, maskPctlOther, pctlOther, wrampProfileOther, candidateStatsOther, figNameOther, cellThirdAdjust, cellBgThird, maskCombinedThird, maskPctlThird, pctlThird, wrampProfileThird, candidateStatsThird, figNameThird)
    % Outputs a figure showing each channel of the image (after 
    % linear contrast ajdustment) amd a figure showing the masks of
    % percentile thresholding, histograms of the intensities in the 
    % cell, and plots the profile.
    
    maskCell = maskCombinedMain > 0;
    if onlyOneChannel == 1
        figure
        subplot
        imshow(cellMainAdjust)
        title(figNameMain, 'FontSize', 8)
    elseif thirdChannel ~= 0
        figure
        subplot(1,3,1)
        imshow(cellMainAdjust)
        title(figNameMain, 'FontSize', 6)

        subplot(1,3,2)
        imshow(cellOtherAdjust)
        title(figNameOther, 'FontSize', 6)

        subplot(1,3,3)
        imshow(cellThirdAdjust)
        title(figNameThird, 'FontSize', 6)
    else
        figure
        subplot(1,2,1)
        imshow(cellMainAdjust)
        title(figNameMain, 'FontSize', 6)

        subplot(1,2,2)
        imshow(cellOtherAdjust)
        title(figNameOther, 'FontSize', 6)
    end
    f1 = gcf;
    
    figure
    if onlyOneChannel == 1
        figure
        subplot(2,3,1)
        imshow(maskCombinedMain, [0 3])
        hold on
        plot(centrCell(1,1), centrCell(1,2), '.','color','blue')
        hold off

        subplot(2,3,2)
        maskPctlConn = bwconncomp(maskPctlMain);
        maskRGB = label2rgb(labelmatrix(maskPctlConn),'hsv','w','shuffle');
        imshow(maskCell)
        hold on
        imshow(maskRGB)
        alpha 0.75
        hold off

        subplot(2,3,3)
        plot(1:length(wrampProfileMain), wrampProfileMain, 'color', 'green');                           
        title('Intensity profile - WRAMP axis', 'FontSize', 4)
        xlim([0 length(wrampProfileMain)])

        subplot(2,3,4)
        histogram(cellBgMain(maskCell == 1), 100)
        hold on
        xline(pctlMain)
        hold off
        title([figNameMain, ' intensity'], 'FontSize', 6)
        set(gca,'FontSize',4)

%         subplot(2,3,5)
%         histogram([candidateStatsMain.sumIntensities], 50)
%         title([figNameMain, ' int. object intensities'], 'FontSize', 5);

%         subplot(2,3,6)
%         boxplot([candidateStatsMain.Area])
%         hold on
%         scatter(ones(1,length([candidateStatsMain.Area])), [candidateStatsMain.Area], '.', 'jitter', 'on')
%         title([figNameMain, ' object areas'], 'FontSize', 5);
%         alpha 0.5;
%         hold off
    else
        if thirdChannel ~= 0
            pos5 = [0.67 0.75, 0.12 0.12];
            subplot('Position', pos5)
            imshow(maskCombinedThird, [0 3])
            hold on
            plot(centrCell(1,1), centrCell(1,2), '.','color','blue')
            hold off

            pos6 = [0.83 0.75 0.12 0.12];
            subplot('Position', pos6)
            maskPctlConn = bwconncomp(maskPctlThird);
            maskRGB = label2rgb(labelmatrix(maskPctlConn),'hsv','w','shuffle');
            imshow(maskCell)
            hold on
            imshow(maskRGB)
            alpha 0.75
            hold off

            pos10 = [0.51 0.45 0.12 0.12];
            subplot('Position', pos10)
            histogram(cellBgThird(maskCell == 1), 100)
            hold on
            xline(pctlThird)
            hold off
            title([figNameThird, ' intensity'], 'FontSize', 6)
            set(gca,'FontSize',4)

%             pos16 = [0.68 0.15 0.1 0.1];
%             subplot('Position', pos16)
%             histogram([candidateStatsThird.sumIntensities], 50)
%             title([figNameThird, ' int. object intensities'], 'FontSize', 5);

%             pos17 = [0.83 0.15 0.1 0.1];
%             subplot('Position', pos17)
%             boxplot([candidateStatsThird.Area])
%             hold on
%             scatter(ones(1,length([candidateStatsThird.Area])), [candidateStatsThird.Area], '.', 'jitter', 'on')
%             title([figNameThird, ' object areas'], 'FontSize', 5);
%             alpha 0.5;
%             hold off

            maskTemp = zeros(size(maskCombinedMain));
            maskTemp(maskCombinedOther == 3) = 1;
            combWrampMask(:,:,1) = maskTemp;
            maskTemp = zeros(size(maskCombinedMain));
            maskTemp(maskCombinedMain == 3) = 1;
            combWrampMask(:,:,2) = maskTemp;
            maskTemp = zeros(size(maskCombinedMain));
            maskTemp(maskCombinedThird == 3) = 1;
            combWrampMask(:,:,3) = maskTemp;
        else
            maskTemp = zeros(size(maskCombinedMain));
            maskTemp(maskCombinedOther == 3) = 1;
            combWrampMask(:,:,1) = maskTemp;
            maskTemp = zeros(size(maskCombinedMain));
            maskTemp(maskCombinedMain == 3) = 1;
            combWrampMask(:,:,2) = maskTemp;
            combWrampMask(:,:,3) = zeros(size(maskCombinedMain));
        end

        pos1 = [0.03 0.75, 0.12 0.12];
        subplot('Position', pos1)
        imshow(maskCombinedMain, [0 3])
        hold on
        plot(centrCell(1,1), centrCell(1,2), '.','color','blue')
        hold off

        pos2 = [0.19 0.75 0.12 0.12];
        subplot('Position', pos2)
        maskPctlConn = bwconncomp(maskPctlMain);
        maskRGB = label2rgb(labelmatrix(maskPctlConn),'hsv','w','shuffle');
        imshow(maskCell)
        hold on
        imshow(maskRGB)
        alpha 0.75
        hold off

        pos3 = [0.35 0.75, 0.12 0.12]; 
        subplot('Position', pos3)
        imshow(maskCombinedOther, [0 3])
        hold on
        plot(centrCell(1,1), centrCell(1,2), '.','color','blue')
        hold off

        pos4 = [0.51 0.75 0.12 0.12];
        subplot('Position', pos4)
        maskPctlConn = bwconncomp(maskPctlOther);
        maskRGB = label2rgb(labelmatrix(maskPctlConn),'hsv','w','shuffle');
        imshow(maskCell)
        hold on
        imshow(maskRGB)
        alpha 0.75
        hold off
        
        pos7 = [0.03 0.45 0.12 0.12];
        subplot('Position', pos7)
        imshow(maskCombinedMain > 0)
        hold on
        imshow(combWrampMask)
        alpha 0.75
        hold off

        pos8 = [0.19 0.45 0.12 0.12];
        subplot('Position', pos8)
        histogram(cellBgMain(maskCell == 1), 100)
        hold on
        xline(pctlMain)
        hold off
        title([figNameMain, ' intensity'], 'FontSize', 6)
        set(gca,'FontSize',4)

        pos9 = [0.35 0.45 0.12 0.12];
        subplot('Position', pos9)
        histogram(cellBgOther(maskCell == 1), 100)
        hold on
        xline(pctlOther)
        hold off
        title([figNameOther, ' intensity'], 'FontSize', 6)
        set(gca,'FontSize',4)

        pos11 = [0.67 0.45 0.12 0.12];
        subplot('Position', pos11)
        hold on
        plot(1:length(wrampProfileMain), wrampProfileMain, 'color', 'green');
        if thirdChannel ~= 0
            plot(1:length(wrampProfileMain), wrampProfileThird, 'color', 'blue');
        end
        plot(1:length(wrampProfileMain), wrampProfileOther, 'color', 'red');
        title('Intensity profile - WRAMP axis', 'FontSize', 4)
        hold off
        xlim([0 length(wrampProfileMain)])
% 
%         pos12 = [0.08 0.15 0.1 0.1];
%         subplot('Position', pos12)
%         histogram([candidateStatsMain.sumIntensities], 50)
%         title([figNameMain, ' int. object intensities'], 'FontSize', 5);

%         pos13 = [0.23 0.15 0.1 0.1];
%         subplot('Position', pos13)
%         boxplot([candidateStatsMain.Area])
%         hold on
%         scatter(ones(1,length([candidateStatsMain.Area])), [candidateStatsMain.Area], '.', 'jitter', 'on')
%         title([figNameMain, ' object areas'], 'FontSize', 5);
%         alpha 0.5;
%         hold off

%         pos14 = [0.38 0.15 0.1 0.1];
%         subplot('Position', pos14)
%         histogram([candidateStatsOther.sumIntensities], 50)
%         title([figNameOther, ' int. object intensities'], 'FontSize', 5);

%         pos15 = [0.53 0.15 0.1 0.1];
%         subplot('Position', pos15)
%         boxplot([candidateStatsOther.Area])
%         hold on
%         scatter(ones(1,length([candidateStatsOther.Area])), [candidateStatsOther.Area], '.', 'jitter', 'on')
%         title([figNameOther, ' object areas'], 'FontSize', 5);
%         alpha 0.5;
%         hold off
    end

    f2 = gcf;
end
    