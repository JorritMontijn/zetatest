function [vecThisDiff,vecThisFrac1,vecThisFrac2,vecThisFracLinear1,vecRefT] = ...
		getTraceOffsetDiff(vecPseudoT1,vecPseudoTrace1,vecPseudoStartT1,vecPseudoT2,vecPseudoTrace2,vecPseudoStartT2,vecRefT,dblUseMaxDur)
	%getTraceOffsetDiff Calculate temporal offset vector difference. Syntax:
	%[vecThisDiff,vecThisFrac1,vecThisFrac2,vecThisFracLinear1,vecRefT] = ...
	%	getTraceOffsetDiff(vecPseudoT1,vecPseudoTrace1,vecPseudoStartT1,vecPseudoT2,vecPseudoTrace2,vecPseudoStartT2,vecRefT,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	%% get data1
	%build interpolated data
	vecEventStartT1 = sort(vecPseudoStartT1);
	[vecRefT1,matTracePerTrial1] = getInterpolatedTimeSeries(vecPseudoT1,vecPseudoTrace1,vecEventStartT1,dblUseMaxDur,vecRefT);
	vecMeanTrace1 = nanmean(matTracePerTrial1,1)';
	vecThisFrac1 = cumsum(vecMeanTrace1) / sum(vecMeanTrace1);
	%vecM = nanmean(matTracePerTrial1,1)';
	%vecMeanTrace1 = cat(1,vecM(1),diff(vecM));
	%vecThisFrac1 = cumsum(vecMeanTrace1) / sum(vecMeanTrace1);
	
	%get linear fractions
	vecThisFracLinear1 = linspace(mean(vecMeanTrace1),sum(vecMeanTrace1),numel(vecMeanTrace1))' / sum(vecMeanTrace1);
	
	%calc deviation
	vecDeviation1 = vecThisFrac1 - vecThisFracLinear1;
	
	%% get data2
	%build interpolated data
	vecEventStartT2 = sort(vecPseudoStartT2);
	[vecRefT2,matTracePerTrial2] = getInterpolatedTimeSeries(vecPseudoT2,vecPseudoTrace2,vecEventStartT2,dblUseMaxDur,vecRefT);
	vecMeanTrace2 = nanmean(matTracePerTrial2,1)';
	vecThisFrac2 = cumsum(vecMeanTrace2) / sum(vecMeanTrace2);
	%vecM = nanmean(matTracePerTrial2,1)';
	%vecMeanTrace2 = cat(1,vecM(1),diff(vecM));
	%vecThisFrac2 = cumsum(vecMeanTrace2) / sum(vecMeanTrace2);
	
	%get linear fractions
	vecThisFracLinear2 = linspace(mean(vecMeanTrace2),sum(vecMeanTrace2),numel(vecMeanTrace2))' / sum(vecMeanTrace2);
	
	%calc deviation
	vecDeviation2 = vecThisFrac2 - vecThisFracLinear2;
	
	%% calc diff
	vecThisDiff = vecDeviation1 - vecDeviation2;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

