function [tempDiff,relSpikeTimesAgg,spikeFrac1,relSpikeTimes1,spikeFrac2,relSpikeTimes2,tempDiffUnSub,fracLinear] = calcTempDiff2(spikesPerTrial1,spikesPerTrial2,useDur,useFastInterp)
% compute temporal offset vectors, syntax:
% [tempDiff,relSpikeTimesAgg,spikeFrac1,relSpikeTimes1,spikeFrac2,relSpikeTimes2,fracDiff,fracLinear] = calcTempDiff2(spikesPerTrial1,spikesPerTrial2,useDur,useFastInterp)
%
% history:
%   v0.9 - 18 February 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% get trial-pooled relative spike times
relSpikeTimes1 = getDistinctSpikes(cell2vec(spikesPerTrial1));
relSpikeTimes2 = getDistinctSpikes(cell2vec(spikesPerTrial2));

%combine, add two artificial spikes
relSpikeTimesAgg = cat(1,relSpikeTimes1,relSpikeTimes2);
if max(relSpikeTimesAgg) < useDur(2)
    relSpikeTimesAgg = [useDur(1);sort(relSpikeTimesAgg(:));useDur(2)];
else
    relSpikeTimesAgg = [useDur(1);sort(relSpikeTimesAgg(:))];
end

%% get temporal difference vector
%cond1 goes to S1_n/T1_n; cond2 goes to S2_n/T2_n
numSp1 = numel(relSpikeTimes1);
numSp2 = numel(relSpikeTimes2);
numEv1 = numel(spikesPerTrial1);
numEv2 = numel(spikesPerTrial2);

%spike fraction #1
uniqueSpikeFracs1 = (1:numSp1)'/numEv1;
if useFastInterp
    spikeFrac1 = lininterp1f([useDur(1);relSpikeTimes1;useDur(2)],[0;uniqueSpikeFracs1;numSp1/numEv1],relSpikeTimesAgg,nan)';
else
    spikeFrac1 = interp1([useDur(1);relSpikeTimes1;useDur(2)],[0;uniqueSpikeFracs1;numSp1/numEv1],relSpikeTimesAgg,'linear');
end
spikeFrac1 = fillnans(spikeFrac1,numSp1,numEv1);

%spike fraction #2
uniqueSpikeFracs2 = (1:numSp2)'/numEv2;
if useFastInterp
    spikeFrac2 = lininterp1f([useDur(1);relSpikeTimes2;useDur(2)],[0;uniqueSpikeFracs2;numSp2/numEv2],relSpikeTimesAgg,nan)';
else
    spikeFrac2 = interp1([useDur(1);relSpikeTimes2;useDur(2)],[0;uniqueSpikeFracs2;numSp2/numEv2],relSpikeTimesAgg,'linear');
end
spikeFrac2 = fillnans(spikeFrac2,numSp2,numEv2);

%take difference
tempDiff = spikeFrac1-spikeFrac2;

%subtract linear
tempDiffUnSub = tempDiff;
fracLinear = tempDiff(1)+(tempDiff(end)-tempDiff(1))*...
    (relSpikeTimesAgg-relSpikeTimesAgg(1))/(relSpikeTimesAgg(end)-relSpikeTimesAgg(1));
tempDiff = tempDiff-fracLinear;

%tempDiff = tempDiff-mean(tempDiff); % mean is subtracted in main latenzy2() function

end