function [peaksRandD,resampD,resampT] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,useDur,resampNum,useParPool,useFastInterp)
%run bootstraps by swapping trials, syntax:
%   [peaksRandD,resampD,resampT] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,useDur,resampNum,useParPool,useFastInterp)
%
% history:
%   v0.6 - 19 February 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025


%% run bootstraps
%swap trials randomly in each resampling
 %#ok<*PFBNS>
resampT = cell(1,resampNum);
resampD = cell(1,resampNum);
peaksRandD = nan(1,resampNum);
spikesAgg = cat(1,spikesPerEvent1,spikesPerEvent2);
idxSpikes0 = cellfun(@numel,spikesAgg)==0;
numEv1 = numel(spikesPerEvent1);
numEv2 = numel(spikesPerEvent2);
numEvTot = numEv1+numEv2;

if useParPool
    parfor resamp=1:resampNum
        %% get random subsample
        % useRand1 = randi(numEvTot,[1,numEv1]);
        % useRand2 = randi(numEvTot,[1,numEv2]);

        %without replacement:
        shuffledIdx = my_randperm(numEvTot);
        useRand1 = shuffledIdx(1:numEv1);
        useRand2 = shuffledIdx(numEv1+1:end);

        spikesPerTrial1_Rand = spikesAgg(useRand1);
        spikesPerTrial2_Rand = spikesAgg(useRand2);
        if all(idxSpikes0(useRand1)) && all(idxSpikes0(useRand2))
            continue;
        end

        %get difference
        [randD,randT] = calcTempDiff2(spikesPerTrial1_Rand,spikesPerTrial2_Rand,useDur,useFastInterp);

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
        %% get random subsample
        % useRand1 = randi(numEvTot,[1,numEv1]);
        % useRand2 = randi(numEvTot,[1,numEv2]);

        %without replacement:
        shuffledIdx = my_randperm(numEvTot);
        useRand1 = shuffledIdx(1:numEv1);
        useRand2 = shuffledIdx(numEv1+1:end);

        spikesPerTrial1_Rand = spikesAgg(useRand1);
        spikesPerTrial2_Rand = spikesAgg(useRand2);
        if all(idxSpikes0(useRand1)) && all(idxSpikes0(useRand2))
            continue;
        end

        %get difference
        [randD,randT] = calcTempDiff2(spikesPerTrial1_Rand,spikesPerTrial2_Rand,useDur,useFastInterp);

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