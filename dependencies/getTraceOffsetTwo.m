function [vecThisDiff,vecThisFrac1,vecThisFrac2] = getTraceOffsetTwo(matTracePerTrial1,matTracePerTrial2)
	%getTraceOffsetTwo Calculate temporal offset vector difference. Syntax:
	%[vecThisDiff,vecThisFrac1,vecThisFrac2] = getTraceOffsetTwo(matTracePerTrial1,matTracePerTrial2)
	%
	%This is a subfunction for zetatstest2().
	
	%%
	%cond1 goes to sum(v_mu1); cond2 goes to sum(v_mu2)
	vecMeanTrace1 = nanmean(matTracePerTrial1,1)';
	vecMeanTrace2 = nanmean(matTracePerTrial2,1)';
	dblMin = min(min(vecMeanTrace1),min(vecMeanTrace2));
	dblMax = max(max(vecMeanTrace1),max(vecMeanTrace2));
	dblRange = (dblMax-dblMin);
	if dblRange == 0
		dblRange = 1;
		warning([mfilename ':ZeroVar'],'Input data has zero variance');
	end
	vecMeanTrace1 = (vecMeanTrace1-dblMin)./dblRange;
	vecMeanTrace2 = (vecMeanTrace2-dblMin)./dblRange;
	
	%get real cumsums
	vecThisFrac1 = cumsum(vecMeanTrace1);
	vecThisFrac2 = cumsum(vecMeanTrace2);
	
	%take difference
	vecDeviation = vecThisFrac1 - vecThisFrac2;
	
	%mean-subtract?
	vecThisDiff = vecDeviation - mean(vecDeviation);
end

