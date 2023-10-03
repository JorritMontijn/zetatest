function [dblZetaP,sZETA] = zetatest2(vecSpikeTimes1,matEventTimes1,vecSpikeTimes2,matEventTimes2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile)
	%zetatest2 Calculates p-value for difference in responsiveness between two neurons
	%syntax: [dblZetaP,sZETA] = zetatest2(vecSpikeTimes1,matEventTimes1,vecSpikeTimes2,matEventTimes2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile)
	%
	%required inputs:
	%	- vecSpikeTimes1 [S x 1]: spike times (in seconds) for neuron 1
	%	- vecEventTimes1 [T x 1]: event on times (s) for neuron 1, or [T x 2] including event off times
	%	- vecSpikeTimes2 [S x 1]: spike times (in seconds) for neuron 2
	%	- vecEventTimes2 [T x 1]: event on times (s) for neuron 2, or [T x 2] including event off times
	%optional inputs:
	%	- boolPairedTest: boolean, use paired or unpaired testing (default=unpaired; boolPairedTest=false). See note below.
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all spikes beyond this duration after event onset
	%								(default: minimum of all event onsets to next event onset)
	%	- intResampNum: integer, number of resamplings (default: 100)
	%	- intPlot: integer, plotting switch (0=none, 1=inst. rate only, 2=traces only, 3=raster plot as well) (default: 0)
	%	- boolDirectQuantile; boolean, switch to use the empirical null-distribution rather than the
	%								Gumbel approximation (default: false) [Note: requires many resamplings!]
	%
	%output:
	%	- dblZetaP; p-value based on Zenith of Event-based Time-locked Anomalies
	%	- sZETA; structure with fields:
	%		- dblZETA; responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecD)
	%		- vecD; temporal deviation vector of data
	%		- cellRandDiff; null-hypothesis temporal deviation vectors of jittered data
	%		- dblD_InvSign; largest peak of inverse sign to ZETA (i.e., -ZETA)
	%		- dblPeakT_InvSign; time corresponding to -ZETA
	%		- intPeakIdx_InvSign; entry corresponding to -ZETA
	%		- dblUseMaxDur; window length used to calculate ZETA
	%
	%*Note on paired versus unpaired testing*
	%The "paired" two-sample zeta-test means the temporal jitters used to estimate the null
	%distribution are identical for the two neurons. This can be useful when two neurons are
	%recorded simultaneously and you wish to test whether they respond differently. With paired
	%testing, any temporally localized noise (e.g., movement or optogenetic artifacts) present
	%in both neurons can cancel each other out and raise the statistical sensitivity. However, this
	%obviously only makes sense if the neurons are recorded simultaneously and when they are equally
	%affected by noise. So to summarize when to use:
	%- Paired: simultaneously recorded & similarly noise-affected neurons from the same experiment
	%- Unpaired: all other cases (e.g., does a neuron respond differently to two stimulus sets?)
	%
	%v0.2 - 11 January 2022
	
	%Version history:
	%0.1 - 7 Dec 2021
	%	Created by Jorrit Montijn
	%0.2 - 11 January 2022
	%	Added paired/unpaired testing, updated syntax & help [by JM]
	%0.3 - 28 Sept 2023
	%	New procedure using trial swapping [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes1 = vecSpikeTimes1(:);
	assert(isnumeric(vecSpikeTimes1),[mfilename ':WrongInputType'], 'Supplied spike time variable 1 is not a numeric vector');
	vecSpikeTimes2 = vecSpikeTimes2(:);
	assert(isnumeric(vecSpikeTimes2),[mfilename ':WrongInputType'], 'Supplied spike time variable 2 is not a numeric vector');
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanZ = nan;
	dblMeanP = nan;
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
		intResampNum = 500;
	end
	
	%get intPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolDirectQuantile
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end

	%% get zeta
	vecEventStarts1 = matEventTimes1(:,1);
	vecEventStarts2 = matEventTimes2(:,1);
	if numel(vecEventStarts1) > 1 && (numel(vecSpikeTimes1)+numel(vecSpikeTimes2)) > 0 && ~isempty(dblUseMaxDur) && dblUseMaxDur>0
		[vecSpikeT,vecRealDiff,vecRealFrac1,vecRealFrac2,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcZetaTwo(vecSpikeTimes1,vecEventStarts1,vecSpikeTimes2,vecEventStarts2,dblUseMaxDur,intResampNum,boolDirectQuantile);
	else
		intZETALoc = nan;
	end
	
	%% build placeholder outputs
	sZETA = [];
	if isnan(intZETALoc)
		dblZetaP = 1;
		dblZETA = 0;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		
		%build placeholder outputs
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = 0;
		sZETA.dblP = 1;
		sZETA.dblPeakT = nan;
		sZETA.intPeakIdx = [];
		if boolStopSupplied
			sZETA.dblMeanD = 0;
			sZETA.dblMeanP = 1;
		end
		sZETA.vecSpikeT = [];
		sZETA.vecD = [];
		sZETA.cellRandDiff = [];
		
		sZETA.dblD_InvSign = 0;
		sZETA.dblPeakT_InvSign = nan;
		sZETA.intPeakIdx_InvSign = [];
		sZETA.dblUseMaxDur = nan;
		return
	end
	
	%% extract real outputs
	%get location
	dblMaxDTime = vecSpikeT(intZETALoc);
	dblD = vecRealDiff(intZETALoc);
	
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblD)*vecRealDiff);
	dblMaxDTimeInvSign = vecSpikeT(intPeakLocInvSign);
	dblD_InvSign = vecRealDiff(intPeakLocInvSign);
	
	%% calculate mean-rate difference with t-test
	if boolStopSupplied && (nargout > 1 || intPlot > 1)
		%neuron 1
		vecRespBinsDur = sort(flat([matEventTimes1(:,1) matEventTimes1(:,2)]));
		vecR = histcounts(vecSpikeTimes1,vecRespBinsDur);
		vecD = diff(vecRespBinsDur)';
		vecMu_Dur1 = vecR(1:2:end)./vecD(1:2:end);
		dblStart1 = min(vecRespBinsDur);
		dblFirstPreDur = dblStart1 - max([0 dblStart1 - median(vecD(2:2:end))]);
		dblR11 = sum(vecSpikeTimes1 > (dblStart1 - dblFirstPreDur) & vecSpikeTimes1 < dblStart1);
		vecMu_Pre1 = [dblR11 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
		
		%neuron 2
		vecRespBinsDur = sort(flat([matEventTimes2(:,1) matEventTimes2(:,2)]));
		vecR = histcounts(vecSpikeTimes2,vecRespBinsDur);
		vecD = diff(vecRespBinsDur)';
		vecMu_Dur2 = vecR(1:2:end)./vecD(1:2:end);
		dblStart2 = min(vecRespBinsDur);
		dblFirstPreDur = dblStart2 - max([0 dblStart2 - median(vecD(2:2:end))]);
		dblR12 = sum(vecSpikeTimes2 > (dblStart2 - dblFirstPreDur) & vecSpikeTimes2 < dblStart2);
		vecMu_Pre2 = [dblR12 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
		
		
		%difference
		vecMu1 = vecMu_Dur1 - vecMu_Pre1;
		vecMu2 = vecMu_Dur2 - vecMu_Pre2;
		
		%get metrics
		[h,dblMeanP,ci,stats]=ttest2(vecMu1,vecMu2);
		dblMeanZ = -norminv(dblMeanP/2);
	end
	
	%% plot
	if intPlot > 1
		%plot maximally 50 traces
		intPlotIters = min([numel(cellRandDiff) 50]);
		
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
		if intPlot > 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes1,vecEventStarts1(:,1),dblUseMaxDur,10000);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Raster plot data 1');
			fixfig;
			grid off;
			
			subplot(2,3,4)
			plotRaster(vecSpikeTimes2,vecEventStarts2(:,1),dblUseMaxDur,10000);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Raster plot data 2');
			fixfig;
			grid off;
		end
		
		subplot(2,3,2)
		plot(vecSpikeT,vecRealFrac1);
		hold on
		plot(vecSpikeT,vecRealFrac2);
		title(sprintf('Real data, data 1 - data 2'));
		xlabel('Time after event (s)');
		ylabel('Fractional position of spike in trial');
		fixfig
		
		subplot(2,3,3)
		cla;
		hold all
		for intIter=1:intPlotIters
			plot(cellRandT{intIter},cellRandDiff{intIter},'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealDiff,'Color',lines(1));
		scatter(dblMaxDTime,vecRealDiff(intZETALoc),'bx');
		scatter(dblMaxDTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time after event (s)');
		ylabel('Deviation difference (\deltas)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), z(Hz)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanZ,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		fixfig
	end
	
	%% build optional output structure
	if nargout > 1
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = dblD;
		sZETA.dblP = dblZetaP;
		sZETA.dblPeakT = dblMaxDTime;
		sZETA.intPeakIdx = intZETALoc;
		if boolStopSupplied
			sZETA.dblMeanZ = dblMeanZ;
			sZETA.dblMeanP = dblMeanP;
			sZETA.vecMu1 = vecMu1;
			sZETA.vecMu2 = vecMu2;
		end
		sZETA.vecSpikeT = vecSpikeT;
		sZETA.vecD = vecRealDiff;
		sZETA.cellRandDiff = cellRandDiff;
		
		sZETA.dblD_InvSign = dblD_InvSign;
		sZETA.dblPeakT_InvSign = dblMaxDTimeInvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblUseMaxDur = dblUseMaxDur;
	end
end

