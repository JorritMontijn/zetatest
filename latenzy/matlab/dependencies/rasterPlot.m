function rasterPlot(spikeTimes, eventTimes, useDur, trialType, plotColor, plotMaxSpikes, addLabels)
% make a raster plot, syntax:
%   rasterPlot(spikeTimes, eventTimes, useDur, trialType, plotColor, plotMaxSpikes, addLabels)
%   inputs:
%   - spikeTimes [S x 1]: spike times (s)*
%   - eventTimes [T x 1]: event (start) times (s)
%   - useDur: scalar or [pre post], time to include after/around event times (s)
%   - trialType: label for trialtype (e.g., orientation), same size as eventTimes
%   - plotColor: colors for plotting [3 x N], where N is unique trialtypes
%   - plotMaxSpikes: max number of spikes to plot (default: inf)
%   - addLabels: boolean, add trialType labels (default: true)
%
% *alternative input available:
%   - spikeTimes should be a cell array, where every cell contains the aligned(!) spikes for a repetition
%   - set eventTimes to []
%   - set useDur to match (or use a smaller window)
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%get useDur
if ~exist('useDur','var') || isempty(useDur)
    sortedEventTimes = sort(eventTimes);
    useDur = min(diff(sortedEventTimes));
end
if isscalar(useDur)
    useDur = sort([0 useDur]);
elseif numel(useDur)~=2
    error([mfilename ':WrongMaxDurInput'],'useDur must be a scalar or a two-element array');
end

if ~exist('trialType', 'var') || isempty(trialType)
    trialType = ones(size(eventTimes));
end
[trialType, uniqueType, ~, ~, ~] = val2idx(trialType);
numTrialType = numel(uniqueType);

if ~exist('plotColor', 'var') || isempty(plotColor) || size(plotColor, 1) ~= numTrialType
    colorsOut = getColors(numTrialType - 1);
    colorsOut = colorsOut(randperm(size(colorsOut, 1)), :);
    plotColor = [0 0 0; colorsOut];
end

if ~exist('plotMaxSpikes', 'var') || isempty(plotMaxSpikes)
    plotMaxSpikes = Inf;
end

if ~exist('addLabels','var') || isempty(addLabels)
    addLabels = true;
end

%% make raster plot
%subselect
if numel(spikeTimes) > plotMaxSpikes
    spikeTimes = spikeTimes(sort(randperm(numel(spikeTimes), plotMaxSpikes)));
end

% cla;
hold on;
offset = 0;
if numTrialType > 1
    for thisTrialType = 1:numel(uniqueType)
        %get event times
        thisTrialStarts = eventTimes(trialType == thisTrialType);

        %get spike times in subset of trials
        [~, spikesPerEvent] = getRelSpikeTimes(spikeTimes, thisTrialStarts, useDur);

        %plot spikes per trial
        for thisTrial = 1:numel(thisTrialStarts)
            theseTimes = spikesPerEvent{thisTrial};
            line([theseTimes(:)'; theseTimes(:)'], [(thisTrial + offset) * ones(1, numel(theseTimes)) - 0.5; (thisTrial + offset) * ones(1, numel(theseTimes)) + 0.5], ...
                'Color', plotColor(thisTrialType, :), 'LineWidth', 1.5);
        end

        %add label for this trial type at half the number of trials
        numTrials = numel(thisTrialStarts);
        yLabelPos = offset + round(numTrials / 2);

        %add label only if there are multiple trial types
        if numTrialType > 1 & addLabels
            text(max(xlim) + 0.05, yLabelPos, sprintf('%d', uniqueType(thisTrialType)), ...
                'Color', plotColor(thisTrialType, :), 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
        offset = offset + numTrials;
    end
else
    [~, spikesPerEvent] = getRelSpikeTimes(spikeTimes, eventTimes, useDur);

    %plot spikes per trial
    for thisTrial = 1:numel(eventTimes)
        theseTimes = spikesPerEvent{thisTrial};
        line([theseTimes(:)'; theseTimes(:)'], [thisTrial * ones(1, numel(theseTimes)) - 0.5; thisTrial * ones(1, numel(theseTimes)) + 0.5], ...
            'Color', plotColor(trialType(thisTrial), :), 'LineWidth', 1.5);
    end
end
hold off;

%set figure properties
ylim([0.5 numel(eventTimes) + 0.5]);
xlim(useDur);
xlabel('Time from event (s)');
ylabel('Trial');

%add labels for each trialType on the right y-axis
if numTrialType > 1 & addLabels

    %add a secondary y-axis for labels
    ax2 = axes('Position', [0.9 0.1 0.05 0.8], 'Color', 'none', 'YAxisLocation', 'right');
    axis(ax2,'off'); % turn off the axis

    for thisTrialType = 1:numTrialType
        numTrials = numel(find(trialType == uniqueType(thisTrialType)));
        yLabelPos = offset + round(numTrials / 2) - numTrials; % Adjust to half
        text(ax2, 0, yLabelPos, sprintf('%d', uniqueType(thisTrialType)), ...
            'Color', plotColor(thisTrialType, :), 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end
end

function colors = getColors(numColors)
% this function generates a list of RGB color codes with the specified
% number of colors using an extended custom palette. The function returns
% an [numColors x 3] matrix of RGB color codes.

%define an extended colorblind-friendly palette
cudPalette = [
    0.0, 0.6, 0.5;
    0.8, 0.6, 0.2;
    0.9, 0.3, 0.4;
    0.5, 0.4, 0.8;
    0.7, 0.7, 0.7;
    0.3, 0.5, 0.8;
    0.8, 0.3, 0.5;
    0.4, 0.7, 0.3;
    0.9, 0.7, 0.1;
    0.5, 0.5, 0.1;
    0.1, 0.5, 0.5;
    0.7, 0.4, 0.2;
    0.6, 0.2, 0.9;
    0.2, 0.6, 0.9;
    0.9, 0.6, 0.4;
    0.4, 0.9, 0.6;
    0.8, 0.8, 0.4;
    0.6, 0.4, 0.1;
    0.3, 0.7, 0.7;
    0.9, 0.4, 0.7;
    0.4, 0.4, 0.4;
    0.6, 0.3, 0.5;
    0.2, 0.8, 0.4;
    0.5, 0.6, 0.7;
    0.8, 0.5, 0.3;
    0.7, 0.2, 0.6;
    0.2, 0.5, 0.8;
    0.9, 0.9, 0.6;
    0.3, 0.3, 0.9;
    0.7, 0.6, 0.3;
    0.6, 0.7, 0.6;
    0.9, 0.5, 0.3;
    0.2, 0.8, 0.6;
    0.5, 0.7, 0.5;
    ];

%ensure there are no more conditions than colors
numPaletteColors = size(cudPalette, 1);
if numColors > numPaletteColors
    error('There are more conditions than colors, please define your own palette');
else
    colors = cudPalette(1:numColors,:);
end
end
