function [peaksRandD,resampD,resampT] = runJitterBootstraps(spikeTimes,eventTimes,useDur,resampNum,jitterSize,useParPool)
% run bootstraps by jittering event times, syntax:
%   [peaksRandD,randD,randT] = runJitterBootstraps(spikeTimes,eventTimes,preEventTime,postEventTime,resampNum,jitterSize,useParPool)
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%ensure correct orientation
spikeTimes = sort(spikeTimes(:));
eventTimes = sort(eventTimes(:));

%% create random bootstraps
resampD = cell(1,resampNum);
resampT = cell(1,resampNum);
peaksRandD = nan(1,resampNum);
eventNum = numel(eventTimes);
fullDuration = useDur(2)-useDur(1);
jitterPerTrial = nan(eventNum,resampNum);

%uniform jitters between jitterSize*[-tau,+tau]
for resamp=1:resampNum
    jitterPerTrial(:,resamp) = jitterSize*fullDuration*((rand(size(eventTimes))-0.5)*2);
end

%% run bootstraps
if useParPool
    parfor resamp=1:resampNum
        randEventT = eventTimes+jitterPerTrial(:,resamp);
        [randD,randT] = calcTempDiff(spikeTimes,randEventT,useDur);

        %get largest deviation
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            maxRandD = minVal;
        else
            maxRandD = maxVal;
        end

        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(maxRandD)
            peaksRandD(resamp) = maxRandD;
        end
    end
else
    for resamp=1:resampNum
        randEventT = eventTimes+jitterPerTrial(:,resamp);
        [randD,randT] = calcTempDiff(spikeTimes,randEventT,useDur);

        %get largest deviation
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            maxRandD = minVal;
        else
            maxRandD = maxVal;
        end
        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(maxRandD)
            peaksRandD(resamp) = maxRandD;
        end
    end
end

end
