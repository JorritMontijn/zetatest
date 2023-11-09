function [vecSpikeT,vecThisDiff,vecThisFrac1,vecSpikeTimes1,vecThisFrac2,vecSpikeTimes2] = ...
		getTempOffsetTwo(cellTimePerSpike1,cellTimePerSpike2,dblUseMaxDur,boolFastInterp)
		%getTempOffsetTwo Calculate temporal offset vectors. Syntax:
	%[vecSpikeT,vecThisDiff,vecThisFrac1,vecSpikeTimes1,vecThisFrac2,vecSpikeTimes2] = ...
	%	getTempOffsetTwo(cellTimePerSpike1,cellTimePerSpike2,dblUseMaxDur,boolFastInterp)
	%This is a subfunction for zetatest2().
	
	%%
	vecSpikeTimes1 = getUniqueSpikes(cell2vec(cellTimePerSpike1));
	vecSpikeTimes2 = getUniqueSpikes(cell2vec(cellTimePerSpike2));
	vecSpikeT = sort(cat(1,vecSpikeTimes1,vecSpikeTimes2));
	vecSpikeT = [0;sort(vecSpikeT(:),'ascend');dblUseMaxDur]; %add start/end
	
	%cond1 goes to S1_n/T1_n; cond2 goes to S2_n/T2_n
	intSp1 = numel(vecSpikeTimes1);
	intSp2 = numel(vecSpikeTimes2);
	intT1 = numel(cellTimePerSpike1);
	intT2 = numel(cellTimePerSpike2);
	
	%spike fraction #1
	vecUniqueSpikeFracs1 = (1:intSp1)'/intT1;
	if boolFastInterp
		vecThisFrac1 = lininterp1f([0;vecSpikeTimes1;dblUseMaxDur],[0;vecUniqueSpikeFracs1;intSp1/intT1],vecSpikeT,nan)';
	else
		vecThisFrac1 = interp1([0;vecSpikeTimes1;dblUseMaxDur],[0;vecUniqueSpikeFracs1;intSp1/intT1],vecSpikeT);
	end
	vecThisFrac1 = fillnans(vecThisFrac1,intSp1,intT1);
	
	%spike fraction #2
	vecUniqueSpikeFracs2 = (1:intSp2)'/intT2;
	if boolFastInterp
		vecThisFrac2 = lininterp1f([0;vecSpikeTimes2;dblUseMaxDur],[0;vecUniqueSpikeFracs2;intSp2/intT2],vecSpikeT,nan)';
	else
		vecThisFrac2 = interp1([0;vecSpikeTimes2;dblUseMaxDur],[0;vecUniqueSpikeFracs2;intSp2/intT2],vecSpikeT);
	end
	vecThisFrac2 = fillnans(vecThisFrac2,intSp2,intT2);
	
	%take difference
	vecDeviation = vecThisFrac1 - vecThisFrac2;
	
	%mean-subtract?
	vecThisDiff = vecDeviation - mean(vecDeviation);
end
