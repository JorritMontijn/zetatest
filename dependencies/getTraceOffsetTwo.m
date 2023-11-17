function [vecThisDiff,vecThisFrac1,vecThisFrac2] = getTraceOffsetTwo(matTracePerTrial1,matTracePerTrial2)
	%getTraceOffsetTwo Calculate temporal offset vector difference. Syntax:
	%[vecThisDiff,vecThisFrac1,vecThisFrac2] = getTraceOffsetTwo(matTracePerTrial1,matTracePerTrial2)
	%
	%This is a subfunction for zetatstest2().
	
	%%
	%cond1 goes to sum(v_mu1); cond2 goes to sum(v_mu2)
	vecMeanTrace1 = nanmean(matTracePerTrial1,1)';
	vecMeanTrace2 = nanmean(matTracePerTrial2,1)';
	
	%get real cumsums
	vecThisFrac1 = cumsum(vecMeanTrace1);
	vecThisFrac2 = cumsum(vecMeanTrace2);
	
	%take difference
	vecDeviation = vecThisFrac1 - vecThisFrac2;
	
	%mean-subtract?
	vecThisDiff = vecDeviation - mean(vecDeviation) + vecDeviation(1);
end

