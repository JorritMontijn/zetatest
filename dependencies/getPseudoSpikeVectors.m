function [vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(vecSpikeTimes,vecEventT,dblWindowDur)
	
	%ensure sorting and alignment
	vecSpikeTimes = sort(vecSpikeTimes(:));
	vecEventT = sort(vecEventT(:));
	
	%% pre-allocate
	intSamples = numel(vecSpikeTimes);
	intTrials = numel(vecEventT);
	dblMedianDur = median(diff(vecSpikeTimes));
	cellPseudoSpikeT = cell(1,intTrials);
	vecPseudoStartT = nan(intTrials,1);
	dblPseudoEventT = 0;
	intLastUsedSample = 0;
	intFirstSample = [];
	%run
	for intTrial=1:intTrials
		dblEventT = vecEventT(intTrial);
		intStartSample = (find(vecSpikeTimes > dblEventT,1));
		intEndSample = (find(vecSpikeTimes > (dblEventT+dblWindowDur),1)-1);
		if intStartSample > intEndSample
			intEndSample = [];
			intStartSample = [];
		end
		if isempty(intEndSample)
			intEndSample = intStartSample;
		end
		vecEligibleSamples = intStartSample:intEndSample;
		indRemSamples = (vecEligibleSamples <= 0) | (vecEligibleSamples > intSamples);
		vecUseSamples = vecEligibleSamples(~indRemSamples);
		
		%check if beginning or end
		if ~isempty(vecUseSamples)
			if intTrial==1
				vecUseSamples = 1:vecUseSamples(end);
			elseif intTrial==intTrials
				vecUseSamples = vecUseSamples(1):intSamples;
			end
		end
		vecAddT = vecSpikeTimes(vecUseSamples);
		indOverlap = (vecUseSamples <= intLastUsedSample);
		
		%get event t
		if intTrial == 1
			dblPseudoEventT = 0;
		else
			if intTrial > 1 && dblWindowDur > (dblEventT - vecEventT(intTrial-1))
				vecUseSamples = vecUseSamples(~indOverlap);
				vecAddT = vecSpikeTimes(vecUseSamples);
				dblPseudoEventT = dblPseudoEventT + dblEventT - vecEventT(intTrial-1);
			else
				dblPseudoEventT = dblPseudoEventT + dblWindowDur;
			end
		end
		
		%% MAKE LOCAL TO EVENT
		if isempty(vecUseSamples)
			vecLocalPseudoT = [];
		else
			intLastUsedSample = vecUseSamples(end);
			vecLocalPseudoT = vecAddT - dblEventT + dblPseudoEventT;
		end
		if isempty(intFirstSample) && ~isempty(vecUseSamples)
			intFirstSample = vecUseSamples(1);
			dblPseudoT0 = dblPseudoEventT;
		end
		
		
		cellPseudoSpikeT{intTrial} = vecLocalPseudoT;
		vecPseudoStartT(intTrial) = dblPseudoEventT;
		
	end
	
	%add beginning
	if intFirstSample > 1
		dblStepBegin = vecSpikeTimes(intFirstSample) - vecSpikeTimes(intFirstSample-1);
		vecSampAddBeginning = 1:(intFirstSample-1);
		cellPseudoSpikeT = cat(2,{vecSpikeTimes(vecSampAddBeginning) - vecSpikeTimes(vecSampAddBeginning(1)) + dblPseudoT0 - dblStepBegin - range(vecSpikeTimes(vecSampAddBeginning))},cellPseudoSpikeT);
	end
	
	%add end
	intTn = numel(vecSpikeTimes);
	if intTn > intLastUsedSample
		vecSampAddEnd = (intLastUsedSample+1):intTn;
		cellPseudoSpikeT = cat(2,cellPseudoSpikeT,{vecSpikeTimes(vecSampAddEnd) - dblEventT + dblPseudoEventT + dblWindowDur});
	end
	
	%recombine into vector
	vecPseudoSpikeTimes = cell2vec(cellPseudoSpikeT);
end

