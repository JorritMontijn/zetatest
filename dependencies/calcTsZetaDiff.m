function [vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,vecRealFracLinear,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTsZetaDiff(vecTraceT1,vecTraceAct1,vecEventStarts1,vecTraceT2,vecTraceAct2,vecEventStarts2,dblSamplingInterval,dblUseMaxDur,intResampNum,boolPairwise,boolDirectQuantile,dblJitterSize)
	%calcTraceZeta Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTsZetaDiff(vecTraceT1,vecTraceAct1,vecEventStarts1,vecTraceT2,vecTraceAct2,vecEventStarts2,dblSamplingInterval,dblUseMaxDur,intResampNum,boolPairwise,boolDirectQuantile,dblJitterSize)
	
	%% check inputs and pre-allocate error output
	vecRefT = [];
	vecRealDiff = [];
	vecRealFrac1 = [];
	vecRealFrac2 = [];
	vecRealFracLinear = [];
	cellRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	
	%% reduce data
	if size(vecEventStarts1,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	%discard leading/lagging data1
	vecEventStarts1 = vecEventStarts1(:,1);
	dblPreUse1 = -dblUseMaxDur*dblJitterSize;
	dblPostUse1 = dblUseMaxDur*(dblJitterSize+1);
	dblStartT1 = min(vecEventStarts1) + dblPreUse1*2;
	dblStopT1 = max(vecEventStarts1) + dblPostUse1*2;
	indRemoveEntries1 = (vecTraceT1 < dblStartT1) | (vecTraceT1 > dblStopT1);
	vecTraceT1(indRemoveEntries1) = [];
	vecTraceAct1(indRemoveEntries1) = [];
	
	%discard leading/lagging data2
	vecEventStarts2 = vecEventStarts2(:,1);
	dblPreUse2 = -dblUseMaxDur*dblJitterSize;
	dblPostUse2 = dblUseMaxDur*(dblJitterSize+1);
	dblStartT2 = min(vecEventStarts2) + dblPreUse2*2;
	dblStopT2 = max(vecEventStarts2) + dblPostUse2*2;
	indRemoveEntries2 = (vecTraceT2 < dblStartT2) | (vecTraceT2 > dblStopT2);
	vecTraceT2(indRemoveEntries2) = [];
	vecTraceAct2(indRemoveEntries2) = [];
	
	%% build pseudo data, stitching stimulus periods
	[vecPseudoT1,vecPseudoTrace1,vecPseudoStartT1] = getPseudoTimeSeries(vecTraceT1,vecTraceAct1,vecEventStarts1,dblUseMaxDur);
	[vecPseudoT2,vecPseudoTrace2,vecPseudoStartT2] = getPseudoTimeSeries(vecTraceT2,vecTraceAct2,vecEventStarts2,dblUseMaxDur);
	vecPseudoTrace1 = vecPseudoTrace1 - min(vecPseudoTrace1(:));
	vecPseudoTrace2 = vecPseudoTrace2 - min(vecPseudoTrace2(:));
	if numel(vecPseudoT1) < 3
		return;
	end
	if numel(vecPseudoT2) < 3
		return;
	end
	
	%% build reference time
	vecRefT1 = getTsRefT(vecPseudoT1,vecPseudoStartT1,dblUseMaxDur);
	vecRefT2 = getTsRefT(vecPseudoT2,vecPseudoStartT2,dblUseMaxDur);
	%set tol
	dblSampInterval = (median(diff(vecPseudoT1)) + median(diff(vecPseudoT2)))/2;
	dblTol = dblSampInterval/100;
	vecRefT = uniquetol(cat(1,vecRefT1(:),vecRefT2(:)),dblTol);
	
	%% get trial responses
	[vecRealDiff,vecRealFrac1,vecRealFrac2,vecRealFracLinear,vecRefT] = ...
		getTraceOffsetDiff(vecPseudoT1,vecPseudoTrace1,vecPseudoStartT1',vecPseudoT2,vecPseudoTrace2,vecPseudoStartT2',vecRefT,dblUseMaxDur);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	%% run bootstraps; try parallel, otherwise run normal loop
	cellRandDiff = cell(1,intResampNum);
	vecMaxRandD = nan(1,intResampNum);
	intTrials1 = numel(vecPseudoStartT1);
	intTrials2 = numel(vecPseudoStartT2);
	%vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %original zeta
	vecJitterPerTrial1 = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials1)'; %new
	matJitterPerTrial1 = nan(intTrials1,intResampNum);
	vecJitterPerTrial2 = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials2)'; %new
	matJitterPerTrial2 = nan(intTrials2,intResampNum);
	for intResampling=1:intResampNum
		matJitterPerTrial1(:,intResampling) = vecJitterPerTrial1(randperm(numel(vecJitterPerTrial1)));
		matJitterPerTrial2(:,intResampling) = vecJitterPerTrial2(randperm(numel(vecJitterPerTrial2)));
	end
	if boolPairwise
		matJitterPerTrial2 = matJitterPerTrial1;
	end
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecPseudoRandEventT1 = vecPseudoStartT1' + matJitterPerTrial1(:,intResampling);
			vecPseudoRandEventT2 = vecPseudoStartT2' + matJitterPerTrial2(:,intResampling);
			
			%% recalc for rand
			vecRandDiff = getTraceOffsetDiff(vecPseudoT1,vecPseudoTrace1,vecPseudoRandEventT1,vecPseudoT2,vecPseudoTrace2,vecPseudoRandEventT2,vecRefT,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	catch
		%normal loop
		for intResampling=1:intResampNum
			%% get random subsample
			vecPseudoRandEventT1 = vecPseudoStartT1' + matJitterPerTrial1(:,intResampling);
			vecPseudoRandEventT2 = vecPseudoStartT2' + matJitterPerTrial2(:,intResampling);
			
			%% recalc for rand
			vecRandDiff = getTraceOffsetDiff(vecPseudoT1,vecPseudoTrace1,vecPseudoRandEventT1,vecPseudoT2,vecPseudoTrace2,vecPseudoRandEventT2,vecRefT,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	end
	
	%% calculate significance
	dblZetaP = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
	%dblZETA = -norminv(((dblZetaP.*2)./(dblZetaP+1))/2);
	%dblZetaP = (1-normcdf(dblZETA))*2;
end

