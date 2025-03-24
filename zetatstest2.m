function [dblZetaP,sZETA] = zetatstest2(vecTime1,vecValue1,matEventTimes1,vecTime2,vecValue2,matEventTimes2,...
		dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblSuperResFactor)
	%zetatstest2 Calculates difference in responsiveness index zeta for two timeseries
	%syntax: [dblZetaP,sZETA] = zetatstest2(vecTime1,vecValue1,matEventTimes1,vecTime2,vecValue2,matEventTimes2,...
	%								dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblSuperResFactor)
	%	input:
	%	- vecTime1 [N x 1]: time (s) corresponding to entries in vecValue1
	%	- vecValue1 [N x 1]: data values (e.g., calcium imaging dF/F0)
	%	- vecEventTimes1 [T x 1]: event on times (s), or [T x 2] including event off times
	%	- vecTime2 [N x 1]: time (s) corresponding to entries in vecValue2
	%	- vecValue2 [N x 1]: data values (e.g., calcium imaging dF/F0)
	%	- vecEventTimes2 [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all values beyond this duration after stimulus onset
	%								[default: minimum of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 250)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=activity heat map as well) (default: 0)
	%	- boolDirectQuantile; boolean, switch to use the empirical null-distribution 
	%									rather than the Gumbel approximation (default: false)
	%	- dblSuperResFactor; scalar, upsampling of data when calculating zeta (default: 100)
	%
	%	output:
	%	- dblZetaP; Zenith of Event-based Time-locked Anomalies: responsiveness p-value
	%	- sZETA; structure with fields:
	%		- dblZetaP: same as first output
	%		- dblZETA; deviation z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblZetaT; time corresponding to ZETA
	%		- intZetaIdx; entry corresponding to ZETA
	%		- dblMeanZ; z-score based on mean-rate stim/base difference (compare with ZETA)
	%		- dblMeanP; p-value based on mean-rate stim/base difference (compare with ZetaP)
	%		- vecMu1; average spiking rate values per event underlying t-test for condition 1
	%		- vecMu2; average spiking rate values per event underlying t-test for condition 2
	%		- vecRefT: timestamps of trace entries (corresponding to vecRealDiff)
	%		- vecRealDiff; real difference in cumulative density of spikes between condition 1&2
	%		- matRandDiff; random differences in cumulative density of spikes
	%
	%v1.0 - rev20231019
	
	%Version history:
	%1.0 - 2023 October 19
	%	Final release candidate [Created by Jorrit Montijn]
	
	%% prep data
	%check inputs
	assert(numel(vecTime1)==numel(vecValue1) && numel(vecTime2)==numel(vecValue2),...
		[mfilename ':InputError'],['Input lengths do not match']);
	assert(min(matEventTimes1(:,1))>=min(vecTime1) && max(matEventTimes1(:,1))<=(max(vecTime1)-dblUseMaxDur) ...
		&& min(matEventTimes2(:,1))>=min(vecTime2) && max(matEventTimes2(:,1))<=(max(vecTime2)-dblUseMaxDur),...
		[mfilename ':InputError'],['Events exist outside of data period']);
	
	%ensure orientation
	vecTime1 = vecTime1(:);
	[vecTime1,vecReorder1] = sort(vecTime1);
	vecValue1 = vecValue1(:);
	vecValue1 = vecValue1(vecReorder1);
	
	vecTime2 = vecTime2(:);
	[vecTime2,vecReorder2] = sort(vecTime2);
	vecValue2 = vecValue2(:);
	vecValue2 = vecValue2(vecReorder2);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanZ = nan;
	dblMeanP = nan;
	matDataPerTrial1 = [];
	matDataPerTrial2 = [];
	vecMu1 = [];
	vecMu2 = [];
	if size(matEventTimes1,2) > 2
		matEventTimes1 = matEventTimes1';
	end
	if size(matEventTimes2,2) > 2
		matEventTimes2 = matEventTimes2';
	end
	if size(matEventTimes1,2) == 2 && size(matEventTimes2,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min([min(diff(matEventTimes1(:,1))) min(diff(matEventTimes2(:,1)))]);
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 250;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolDirectQuantile
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%get dblJitterSize
	if ~exist('dblSuperResFactor','var') || isempty(dblSuperResFactor)
		dblSuperResFactor = 100; %original:100
	end
	
	%% get ts-zeta diff
	vecEventStarts1 = matEventTimes1(:,1);
	vecEventStarts2 = matEventTimes2(:,1);
	if numel(vecEventStarts1) > 1 && numel(vecTime1) > 1 && numel(vecEventStarts2) > 1 && numel(vecTime2) > 1 && ~isempty(dblUseMaxDur) && dblUseMaxDur>0
		[vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZetaIdx,matDataPerTrial1,matDataPerTrial2] = ...
			calcTsZetaTwo(vecTime1,vecValue1,vecEventStarts1,vecTime2,vecValue2,vecEventStarts2,dblSuperResFactor,dblUseMaxDur,intResampNum,boolDirectQuantile);
	else
		intZetaIdx = nan;
	end
	
	%get location
	dblZetaT = vecRefT(intZetaIdx);
	dblD = vecRealDiff(intZetaIdx);
	
	%% build placeholder outputs
	sZETA = [];
	if isnan(intZetaIdx)
		dblZetaP = 1;
		dblZETA = 0;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		
		%build placeholder outputs
		sZETA = struct;
		sZETA.dblZetaP = dblZetaP;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = [];
		sZETA.dblZetaT = [];
		sZETA.intZetaIdx = [];
		if boolStopSupplied
			sZETA.dblMeanZ = [];
			sZETA.dblMeanP = [];
			sZETA.vecMu1 = [];
			sZETA.vecMu2 = [];
		end
		sZETA.vecRefT = [];
		sZETA.vecRealDiff = [];
		sZETA.matRandDiff = [];
		sZETA.matDataPerTrial1 = [];
		sZETA.matDataPerTrial2 = [];
		return
	end
	
	%% calculate mean-rate difference with t-test
	if boolStopSupplied
		for intTrace=1:2
			if intTrace == 1
				%pre-allocate
				vecEventStarts = matEventTimes1(:,1);
				vecEventStops = matEventTimes1(:,2);
				vecThisTraceT = vecTime1;
				vecThisTraceAct = vecValue1;
			else
				%pre-allocate
				vecEventStarts = matEventTimes2(:,1);
				vecEventStops = matEventTimes2(:,2);
				vecThisTraceT = vecTime2;
				vecThisTraceAct = vecValue2;
			end
			intTimeNum = numel(vecThisTraceT);
			intMaxRep = numel(vecEventStarts);
			vecStimAct = nan(1,intMaxRep);
			
			%go through trials to build spike time vector
			for intEvent=1:intMaxRep
				%% get original times
				dblStimStartT = vecEventStarts(intEvent);
				dblStimStopT = vecEventStops(intEvent);
				dblBaseStopT = dblStimStartT + dblUseMaxDur;
				
				intStartT = max([1 find(vecThisTraceT > dblStimStartT,1) - 1]);
				intStopT = min([intTimeNum find(vecThisTraceT > dblStimStopT,1) + 1]);
				intEndT = min([intTimeNum find(vecThisTraceT > dblBaseStopT,1) + 1]);
				vecSelectFramesBase = (intStopT+1):intEndT;
				vecSelectFramesStim = intStartT:intStopT;
				
				%% get data
				vecUseStimTrace = vecThisTraceAct(vecSelectFramesStim);
				
				%% get activity
				vecStimAct(intEvent) = mean(vecUseStimTrace);
			end
			
			if intTrace == 1
				vecMu_Dur1 = vecStimAct;
			else
				vecMu_Dur2 = vecStimAct;
			end
		end
		
		%difference
		vecMu1 = vecMu_Dur1;
		vecMu2 = vecMu_Dur2;
		
		%get metrics
		[h,dblMeanP,ci,stats]=ttest2(vecMu1,vecMu2);
		dblMeanZ = -norminv(dblMeanP/2);
	end
	
	%% plot
	if intPlot
		%plot maximally 100 traces
		intPlotIters = min([size(matRandDiff,2) 100]);
		
		%maximize figure
		figure;
		drawnow;
		try
			try
				%try new method
				h = handle(gcf);
				h.WindowState = 'maximized';
			catch
				%try old method with javaframe (deprecated as of R2021)
				sWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				drawnow;
				warning(sWarn);
			end
		catch
		end
		if intPlot > 1
			subplot(2,3,1)
			imagesc(vecRefT,1:size(matDataPerTrial1,1),matDataPerTrial1,[0 1]);
			colormap(hot);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Rescaled data, condition 1');
			colorbar;
			grid off;
			
			subplot(2,3,4)
			imagesc(vecRefT,1:size(matDataPerTrial2,1),matDataPerTrial2,[0 1]);
			colormap(hot);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Rescaled data, condition 2');
			colorbar;
			grid off;
		end
		
		subplot(2,3,2)
		plot(vecRefT,vecRealFrac1);
		hold on
		plot(vecRefT,vecRealFrac2);
		title(sprintf('Real data, data 1 - data 2'));
		xlabel('Time after event (s)');
		ylabel('Cumululative sum in trial');
		xlim([0 dblUseMaxDur]);
		
		subplot(2,3,3)
		cla;
		hold all
		for intIter=1:intPlotIters
			plot(vecRefT,matRandDiff(intIter,:),'Color',[0.5 0.5 0.5]);
		end
		plot(vecRefT,vecRealDiff,'Color',lines(1));
		hold off
		xlabel('Time after event (s)');
		ylabel('Deviation difference (\deltas)');
		if boolStopSupplied
			title(sprintf('T-ZETA2=%.3f (p=%.3f), z(mean)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanZ,dblMeanP));
		else
			title(sprintf('T-ZETA2=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		xlim([0 dblUseMaxDur]);
		fixfig
	end
	
	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZetaP = dblZetaP;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = dblD;
		sZETA.dblZetaT = dblZetaT;
		sZETA.intZetaIdx = intZetaIdx;
		if boolStopSupplied
			sZETA.dblMeanZ = dblMeanZ;
			sZETA.dblMeanP = dblMeanP;
			sZETA.vecMu1 = vecMu1;
			sZETA.vecMu2 = vecMu2;
		end
		sZETA.vecRefT = vecRefT;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
		sZETA.matDataPerTrial1 = matDataPerTrial1;
		sZETA.matDataPerTrial2 = matDataPerTrial2;
	end
end
