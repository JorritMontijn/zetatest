function figHandles = makeLatenzyFigs(sLatenzy,spikeTimes,eventTimes,useDur,makePlots)
% make figures for latenzy method
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%#ok<*AGROW>
latency = sLatenzy.latency;
peakTimes = sLatenzy.peakTimes;
peakVals = sLatenzy.peakVals;
realFrac = sLatenzy.realFrac;
fracLin = sLatenzy.fracLin;
realDiff = sLatenzy.realDiff;
realTime = sLatenzy.realTime;
meanRealDiff = sLatenzy.meanRealDiff;
randDiff = sLatenzy.randDiff;
randTime = sLatenzy.randTime;
meanRandDiff = sLatenzy.meanRandDiff;
pValPeaks = sLatenzy.pValsPeak;
% peakZ = sLatenzy.peakZ;
latenzyIdx = sLatenzy.latenzyIdx;

numIters = numel(peakTimes);
useColors = getAbyss(numIters);
lineWidth = 1.5;
markerSize = 60; %50;

%% PLOT
figure;
figHandles = nan(1,4);
if makePlots==1
    %raster plot
    figHandles(1) = subplot(2,2,1);
    rasterPlot(spikeTimes,eventTimes,useDur);
    yLim = ylim;
    hold on
    plot([latency latency], yLim,'color',[0.8627 0.0784 0.2353],'LineStyle','--','LineWidth',lineWidth);
    set(gca,'box','off','TickDir','out');
    ylabel('Trial');
    xlabel('Time from event (s)');
    title('Aligned spikes');
end

%plot cumulative spikes over time
figHandles(2) = subplot(2,2,2); hold on
for iter = 1:numIters
    plot(realTime{iter},fracLin{iter},'color',[0.5 0.5 0.5],'LineWidth',lineWidth);
    p(iter) =  plot(realTime{iter},realFrac{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
    labels{iter} = sprintf('%d', iter);
end
set(gca,'box','off','TickDir','out');
xlim(useDur);
xlabel('Time from event (s)');
ylabel('Fractional spike position');
title('Cumulative spikes');
lgd = legend(p, labels, 'Location', 'southeast','Box','off');
title(lgd, 'Iteration');

%plot offset from linear baseline
figHandles(3) = subplot(2,2,3); hold on
plot(useDur,[0 0],'color',[0.5 0.5 0.5],'LineWidth',lineWidth,'LineStyle','--');
for iter = 1:numIters
    plot(realTime{iter},realDiff{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
end
scatter(peakTimes(~latenzyIdx),peakVals(~latenzyIdx),markerSize,'x','MarkerEdgeColor',[0 0 0],'LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),peakVals(latenzyIdx),markerSize,'x','MarkerEdgeColor',[0.8627 0.0784 0.2353],'LineWidth',lineWidth);
set(gca,'box','off','TickDir','out');
xlim(useDur);
xlabel('Time from event (s)');
ylabel('Deviation (Δfraction)');
title('Deviation from uniform');

%for iteration w/latency, plot real + jitters
figHandles(4) = subplot(2,2,4); hold on
for thisShuffle=1:length(randDiff)
    plot(randTime{thisShuffle,latenzyIdx},(randDiff{thisShuffle,latenzyIdx}-meanRandDiff(thisShuffle,latenzyIdx)) ,'Color',[0.5 0.5 0.5 0.5],'LineWidth',lineWidth);
end
plot(realTime{latenzyIdx},(realDiff{latenzyIdx}-meanRealDiff(latenzyIdx)),'color',useColors(latenzyIdx,:),'LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),(peakVals(latenzyIdx)-meanRealDiff(latenzyIdx)),markerSize,'x','MarkerEdgeColor',[0.8627 0.0784 0.2353],'LineWidth',lineWidth);
if latenzyIdx(1),xlim(useDur);
else, xlim([useDur(1)  peakTimes(find(latenzyIdx)-1)]);end
set(gca,'box','off','TickDir','out');
xlabel('Time from event (s)');
ylabel('Deviation (Δfraction)');
title(sprintf('Real + jittered data (p=%.4f)',pValPeaks(latenzyIdx)));

%add title
sgtitle(sprintf('latenZy estimate = %.4fs', latency), 'FontWeight', 'bold');
end