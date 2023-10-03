function [vecSpikeT,vecThisDiff,vecThisFrac1,vecSpikeTimes1,vecThisFrac2,vecSpikeTimes2] = ...
		getTempOffsetTwo(cellTimePerSpike1,cellTimePerSpike2,dblUseMaxDur)
		%getTempOffsetTwo Calculate temporal offset vectors. Syntax:
	%[vecSpikeT,vecThisDiff,vecThisFrac1,vecSpikeTimes1,vecThisFrac2,vecSpikeTimes2] = ...
	%	getTempOffsetTwo(cellTimePerSpike1,cellTimePerSpike2,dblUseMaxDur)
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
	dblR1 = intSp1/intT1;
	dblR2 = intSp2/intT2;
	
	%spike fraction #1
	vecUniqueSpikeFracs1 = (1:intSp1)'/intT1;
	vecThisFrac1 = interp1([0;vecSpikeTimes1;dblUseMaxDur],[0;vecUniqueSpikeFracs1;intSp1/intT1],vecSpikeT);
	vecThisFrac1 = fillnans(vecThisFrac1,intSp1,intT1);
	%get linear fractions
	vecThisFracLin1 = (vecSpikeT./dblUseMaxDur)*dblR1;
	vecThisDev1 = vecThisFrac1-vecThisFracLin1;
	
	%spike fraction #2
	vecUniqueSpikeFracs2 = (1:intSp2)'/intT2;
	vecThisFrac2 = interp1([0;vecSpikeTimes2;dblUseMaxDur],[0;vecUniqueSpikeFracs2;intSp2/intT2],vecSpikeT);
	vecThisFrac2 = fillnans(vecThisFrac2,intSp2,intT2);
	%get linear fractions
	vecThisFracLin2 = (vecSpikeT./dblUseMaxDur)*dblR2;
	vecThisDev2 = vecThisFrac2-vecThisFracLin2;
	
	%take difference
	vecDeviation = vecThisDev1*dblR1 - vecThisDev2*dblR2;
	
	%mean-subtract?
	vecThisDiff = vecDeviation - mean(vecDeviation);
end
