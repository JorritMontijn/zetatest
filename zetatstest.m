function [dblZetaP,sZETA] = zetatstest(vecTime,vecValue,matEventTimes,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize,boolUseSuperResolution)
	%zetatstest Calculates responsiveness index zeta for timeseries data
	%syntax: [dblZetaP,sZETA] = zetatstest(vecTime,vecValue,vecEventTimes,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize,boolUseSuperResolution)
	%	input:
	%	- vecTime [N x 1]: time (s) corresponding to entries in vecValue
	%	- vecValue [N x 1]: data values (e.g., calcium imaging dF/F0)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all values beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 100)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=activity heat map as well) (default: 0)
	%	- boolDirectQuantile; boolean, switch to use the empirical
	%							null-distribution rather than the Gumbel approximation (default: false)
	%	- dblJitterSize: scalar, sets the temporal jitter window relative to dblUseMaxDur (default: 2)
	%	- boolUseSuperResolution: boolean, use super-resolution procedure (default: 1)
	%
	%	output:
	%	- dblZetaP; Zenith of Event-based Time-locked Anomalies: responsiveness z-score (i.e., >2 is significant)
	%	- sZETA; structure with fields:
	%		- dblZETA; FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecTraceT: timestamps of trace entries (corresponding to vecZ)
	%		- vecD; temporal deviation vector of data
	%		- matRandD; baseline temporal deviation matrix of jittered data
	%
	%Version history:
	%1.0 - 2021 October 29
	%	Created by Jorrit Montijn
	
	%% prep data
	%ensure orientation
	vecTime = vecTime(:);
	[vecTime,vecReorder] = sort(vecTime);
	vecValue = vecValue(:);
	vecValue = vecValue(vecReorder);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanZ = nan;
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	if size(matEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(matEventTimes(:,1)));
	end
	if numel(dblUseMaxDur)>1
		dblUseMaxDurTtest = dblUseMaxDur(2);
		dblUseMaxDur = dblUseMaxDur(1);
	end
		
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 250;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolVerbose
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%get dblJitterSize
	if ~exist('dblJitterSize','var') || isempty(dblJitterSize)
		dblJitterSize = 2; %original:1
	end
	
	%get boolUseSuperResolution
	if ~exist('boolUseSuperResolution','var') || isempty(boolUseSuperResolution)
		boolUseSuperResolution = 1;
	end
	
	%sampling interval
	dblSamplingInterval = median(diff(vecTime));
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	
	%% gettacezeta
	if boolUseSuperResolution == 0
		[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcTsZeta(vecTime,vecValue,vecEventStarts,dblSamplingInterval,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize);
	elseif boolUseSuperResolution == 1
		[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcTsZetaSR_interp(vecTime,vecValue,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize);
	elseif boolUseSuperResolution == 2
		[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcTsZetaSR(vecTime,vecValue,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize);
	end
	%get location
	dblMaxDTime = vecRefT(intZETALoc);
	dblD = vecRealDiff(intZETALoc);
	
	%% calculate mean-rate difference
	intMaxRep = size(vecEventStarts,1);
	vecStimAct = zeros(intMaxRep,1);
	vecBaseAct = zeros(intMaxRep,1);
	if boolStopSupplied
		%pre-allocate
		vecEventStops = matEventTimes(:,2);
		intTimeNum = numel(vecTime);
		if ~exist('dblUseMaxDurTtest','var') || ~isempty(dblUseMaxDurTtest)
			dblUseMaxDurTtest = dblUseMaxDur;
		end
		
		%go through trials to build spike time vector
		for intEvent=1:intMaxRep
			%% get original times
			dblStimStartT = vecEventStarts(intEvent);
			dblStimStopT = vecEventStops(intEvent);
			dblBaseStopT = dblStimStartT + dblUseMaxDurTtest;
			
			intStartT = max([1 find(vecTime > dblStimStartT,1) - 1]);
			intStopT = min([intTimeNum find(vecTime > dblStimStopT,1) + 1]);
			intEndT = min([intTimeNum find(vecTime > dblBaseStopT,1) + 1]);
			vecSelectFramesBase = (intStopT+1):intEndT;
			vecSelectFramesStim = intStartT:intStopT;
			
			%% get data
			vecUseBaseTrace = vecValue(vecSelectFramesBase);
			vecUseStimTrace = vecValue(vecSelectFramesStim);
			
			%% get activity
			vecBaseAct(intEvent) = mean(vecUseBaseTrace);
			vecStimAct(intEvent) = mean(vecUseStimTrace);
		end
		
		%get metrics
		indUseTrials = ~isnan(vecStimAct) & ~isnan(vecBaseAct);
		vecMu_Dur = vecStimAct(indUseTrials);
		vecMu_Pre = vecBaseAct(indUseTrials);
		[h,dblMeanP]=ttest(vecMu_Dur,vecMu_Pre);
		dblMeanZ = -norminv(dblMeanP/2);
	end
	
	%% plot
	if intPlot
		%plot maximally 100 traces
		intPlotIters = min([numel(cellRandDiff) 100]);
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot > 1
			[vecRefT2,matTracePerTrial] = getTraceInTrial(vecTime,vecValue,vecEventStarts,dblSamplingInterval,dblUseMaxDur);
			subplot(2,3,1)
			imagesc(vecRefT2,1:size(matTracePerTrial,1),matTracePerTrial);
			colormap(hot);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Z-scored activation');
			fixfig;
			grid off;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		sOpt.vecWindow = [0 dblUseMaxDur];
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecTime,vecValue,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		%ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean value over trials'));
		xlabel('Time after event (s)');
		ylabel('Trace value');
		fixfig
		
		subplot(2,3,3)
		plot(vecRefT,vecRealFrac)
		hold on
		plot(vecRefT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
		hold off
		title(sprintf('Real data'));
		xlabel('Time after event (s)');
		ylabel('Fractional position of value in trial');
		fixfig
		
		subplot(2,3,4)
		cla;
		hold all
		for intIter=1:intPlotIters
			plot(cellRandT{intIter},cellRandDiff{intIter},'Color',[0.5 0.5 0.5]);
		end
		plot(vecRefT,vecRealDiff,'Color',lines(1));
		hold off
		xlabel('Time after event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), z(mean)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanZ,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		fixfig
		
		if intPlot > 1 && boolUseSuperResolution == 1
			%set tol
			dblSampInterval = median(diff(vecTime));
			dblTol = dblSampInterval/100;
			vecRef2T = uniquetol(vecRefT,dblTol);
			
			%build interpolated data
			[vecRef3T,matTracePerTrialSR] = getInterpolatedTimeSeries(vecTime,vecValue,vecEventStarts(:,1),dblUseMaxDur,vecRef2T);
			indRemPoints = vecRef3T<0 | vecRef3T>dblUseMaxDur;
			vecRef3T(indRemPoints) = [];
			matTracePerTrialSR(:,indRemPoints)=[];
			vecMeanTrace = nanmean(matTracePerTrialSR,1)';
			
			subplot(2,3,5)
			imagesc(vecRef3T,1:size(matTracePerTrialSR,1),matTracePerTrialSR);
			colormap(hot);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			fixfig;
			grid off;
			
			subplot(2,3,6)
			plot(vecRef3T,vecMeanTrace);
			xlabel('Time after event (s)');
			ylabel('Data value');
			fixfig;
			grid off;
		end
	end
	
	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = dblD;
		sZETA.dblP = dblZetaP;
		sZETA.dblPeakT = dblMaxDTime;
		sZETA.intPeakIdx = intZETALoc;
		if boolStopSupplied
			sZETA.dblMeanZ = dblMeanZ;
			sZETA.dblMeanP = dblMeanP;
			sZETA.vecMu_Dur = vecMu_Dur;
			sZETA.vecMu_Pre = vecMu_Pre;
		end
		sZETA.vecTime = vecTime;
		sZETA.vecD = vecRealDiff;
		sZETA.cellRandT = cellRandT;
		sZETA.cellRandDiff = cellRandDiff;
		sZETA.vecMeanBase = vecStimAct;
		sZETA.vecMeanStim = vecBaseAct;
	end
end
