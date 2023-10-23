function [dblZetaP,sZETA,sRate,sLatencies] = zetatest(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampNum,intPlot,vecRestrictRange,boolDirectQuantile,dblJitterSize,boolStitch)
	%zetatest Calculates neuronal responsiveness index zeta
	%syntax: [dblZetaP,sZETA,sRate,vecLatencies] = zetatest(vecSpikeTimes,vecEventTimes,dblUseMaxDur,intResampNum,intPlot,vecRestrictRange,boolDirectQuantile,dblJitterSize,boolStitch)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (in seconds)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times to calculate mean-rate difference
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all spikes beyond this duration after event onset
	%							[default: minimum of all event onsets to next event onset]
	%	- intResampNum: integer, number of resamplings (default: 100)
	%							[Note: if your p-value is close to significance, you should increase this number to enhance the precision]
	%	- intPlot: integer, plotting switch (0=none, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot) (default: 0)
	%	- vecRestrictRange: temporal range within which to restrict onset/peak latencies (default: [-inf inf])
	%	- boolDirectQuantile; boolean, switch to use the empirical null-distribution rather than the
	%							Gumbel approximation (default: false) [Note: requires many resamplings!]
	%	- dblJitterSize; scalar, sets the temporal jitter window relative to dblUseMaxDur (default: 2)
	%	- boolStitch; boolean, use data-stitching to ensure continuous time (default: true)
	%
	%	output:
	%	- dblZetaP; p-value based on Zenith of Event-based Time-locked Anomalies
	%	- sZETA; structure with fields:
	%		- dblZETA; responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblZetaT; time corresponding to ZETA
	%		- intZetaIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecD)
	%		- vecD; temporal deviation vector of data
	%		- matRandD; baseline temporal deviation matrix of jittered data
	%		- dblD_InvSign; largest peak of inverse sign to ZETA (i.e., -ZETA)
	%		- dblZetaT_InvSign; time corresponding to -ZETA
	%		- intZetaIdx_InvSign; entry corresponding to -ZETA
	%		- dblUseMaxDur; window length used to calculate ZETA
	%		- vecLatencies; different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:
	%			1) Latency of ZETA
	%			2) Latency of largest z-score with inverse sign to ZETA
	%			3) Peak time of instantaneous firing rate
	%			4) Onset time of above peak, defined as the first crossing of peak half-height
	%		- vecLatencyVals; data values corresponding to the above peaks
	%	- sRate; structure with fields:
	%		- vecRate; instantaneous spiking rates (like a PSTH)
	%		- vecT; time-points corresponding to vecRate (same as sZETA.vecSpikeT)
	%		- vecM; Mean of multi-scale derivatives
	%		- vecScale; timescales used to calculate derivatives
	%		- matMSD; multi-scale derivatives matrix
	%		- vecV; values on which vecRate is calculated (same as sZETA.vecZ)
	%		Data on the peak: (only if intLatencyPeaks > 0)
	%		- dblPeakTime; time of peak (in seconds)
	%		- dblPeakWidth; duration of peak (in seconds)
	%		- vecPeakStartStop; start and stop time of peak (in seconds)
	%		- intPeakLoc; spike index of peak (corresponding to sZETA.vecSpikeT)
	%		- vecPeakStartStopIdx; spike indices of peak start/stop (corresponding to sZETA.vecSpikeT)
	%		Additionally, it will return peak onset latency (first crossing of peak half-height) using getOnset.m:
	%		- dblOnset: latency for peak onset
	%	- sLatencies; structure containing latencies (copy of sZETA.vecLatencies)
	%		- Onset
	%		- Peak
	%		- ZETA
	%		- ZETA_InvSign
	%
	%Note: zetatest will use parallel computing if you have an active worker pool; if not, it will
	%not start a parallel pool itself.
	%
	%v3.8 - rev20231023
	
	%Version history:
	%0.9 - 27 June 2019
	%	Created by Jorrit Montijn
	%1.0 - 24 September 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - 27 January 2020
	%	New peak detection procedure using multi-scale derivatives [by JM]
	%2.1 - 5 February 2020
	%	Minor changes and bug fixes [by JM]
	%2.2 - 11 February 2020
	%	Peak width, analytical ZETA correction [by JM]
	%2.3 - 26 February 2020
	%	MSD-based instantaneous spiking rates, onset latency [by JM]
	%2.4 - 11 March 2020
	%	Statistical significance using Gumbel approximation [by	JM]
	%2.5 - 27 May 2020
	%	Standardized syntax and variable names [by JM]
	%2.6 - 27 Nov 2020
	%	Improved computation time; now uses parallel bootstrapping [by JM]
	%	In case of only requesting dblZetaP, computation is now up to 10x faster
	%2.7 - 21 Jan 2021
	%	Improved computation time, using calcZeta/getSpikeT subfunctions [by JM]
	%2.8 - 23 Sept 2021
	%	Added switch to use empirical null distribution for significance calculation [by JM]
	%2.9 - 29 Oct 2021
	%	Added option to change the jitter size [by JM]
	%2.10 - 3 Nov 2021
	%	Fixed figure maximization in new matlab versions that deprecated javaframe functionality [by JM]
	%3.0 - 14 Dec 2021
	%	Added stitching, changed name [by JM]
	%3.1 - 11 Jan 2022
	%	Updated syntax [by JM]
	%3.2 - 2 Mar 2022
	%	Fixed stitching bug for low spiking rates & variable ITIs; and added stitching switch [by JM]
	%3.2.1 - 5 Dec 2022
	%	Bug fix when matEventTimes is empty [by JM]
	%3.3 - 26 May 2023
	%	Faster computation time for IFR calculation - parfor enabled by default when parpool is active [by JM]
	%3.4 - 30 May 2023
	%	Changed parallel pool behaviour for zetatest as well: parfor enabled by default when parpool is active [by JM]
	%	Small changes:  getTempOffsetOne.m now adds minimal jitter to duplicate spikes rather than removing them,
	%					getMultiScaleDeriv.m now discards the artificial begin/end points
	%					zetatest now calculates the IFR when the 3rd output (sRate) is requested
	%3.5 - 6 June 2023
	%	Fixed IFR bug (vecTime/vecRate mismatch) introduced in 3.4, and fixed temporal asymmetry in IFR calculation [by JM]
	%3.5.1 - 24 July 2023
	%	Added isaxes() dependency file, and fixed parallel processing bug for non-windows systems [by JM]
	%3.5.2 - 21 August 2023
	%   Made some minor changes to produce deterministic output identical to the python version
	%   Fixed a bug in getPseudoSpikeVectors that could discard some spikes in the final trial [by JM]
	%3.6 - 23 August 2023
	%   Fixed plotting of onsets which apparently broke some time in the past... [by JM]
	%3.6.1 - 25 August 2023
	%   Added time-sorting step to vecSpikeTimes in case spikes are supplied in random order [by JM]
	%3.7 - 19 Oct 2023
	%   Updated documentation [by JM]
	%3.8 - 23 Oct 2023
	%   Moved plotting to separate subfunction [by JM]
	
	%% prep data
	%ensure orientation
	assert(isnumeric(vecSpikeTimes),[mfilename ':WrongInputType'], 'Supplied spike time variable is not a numeric vector');
	vecSpikeTimes = sort(vecSpikeTimes(:));
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanZ = nan;
	dblMeanP = nan;
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	if size(matEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min(diff(matEventTimes(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 100;%100;
	end
	
	%get intPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get vecRestrictRange
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [-inf inf];
	end
	
	%get boolDirectQuantile
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%get dblJitterSize
	if ~exist('dblJitterSize','var') || isempty(dblJitterSize)
		dblJitterSize = 2;
	end
	assert(dblJitterSize>0,[mfilename ':WrongJitterInput'], sprintf('dblJitterSize must be a positive scalar, you requested %.3f',dblJitterSize));
	
	%get boolStitch
	if ~exist('boolStitch','var') || isempty(boolStitch)
		boolStitch = true;
	end
	
	%% get zeta
	if numel(matEventTimes) > 1 && numel(vecSpikeTimes) > 1 && ~isempty(dblUseMaxDur) && dblUseMaxDur>0
		vecEventStarts = matEventTimes(:,1);
		boolUseParallel = [];
		[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcZetaOne(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolStitch,boolUseParallel);
	else
		intZETALoc = nan;
		vecEventStarts = [];
	end
	
	%% build placeholder outputs
	sZETA = [];
	sRate = [];
	sLatencies = [];
	vecLatencies = nan(1,4);
	vecLatencyVals = nan(1,4);
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
		sZETA.vecRealFrac = [];
		sZETA.vecRealFracLinear = [];
		sZETA.vecD = [];
		sZETA.matRandD = [];
		
		sZETA.dblD_InvSign = 0;
		sZETA.dblPeakT_InvSign = nan;
		sZETA.intPeakIdx_InvSign = [];
		sZETA.dblUseMaxDur = nan;
		sZETA.vecLatencies = vecLatencies;
		sZETA.vecLatencyVals = vecLatencyVals;
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
		vecRespBinsDur = sort(flat([matEventTimes(:,1) matEventTimes(:,2)]));
		vecR = histcounts(vecSpikeTimes,vecRespBinsDur);
		vecD = diff(vecRespBinsDur)';
		vecMu_Dur = vecR(1:2:end)./vecD(1:2:end);
		dblStart1 = min(vecRespBinsDur);
		dblFirstPreDur = dblStart1 - max([0 dblStart1 - median(vecD(2:2:end))]) + eps;
		dblR1 = sum(vecSpikeTimes > (dblStart1 - dblFirstPreDur) & vecSpikeTimes < dblStart1);
		vecMu_Pre = [dblR1 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
		
		%get metrics
		[h,dblMeanP,ci,stats]=ttest(vecMu_Dur,vecMu_Pre);
		dblMeanZ = -norminv(dblMeanP/2);
	end
	
	%% calculate instantaneous firing rates
	if nargout > 2 || intPlot > 0
		%get average of multi-scale derivatives, and rescaled to instantaneous spiking rate
		dblMeanRate = (numel(vecSpikeT)/(dblUseMaxDur*numel(vecEventStarts)));
		[vecRate,sRate] = getMultiScaleDeriv(vecSpikeT,vecRealDiff,[],[],[],dblMeanRate,dblUseMaxDur);
	end
	
	%% calculate IFR statistics
	if ~isempty(sRate)
		%get IFR peak
		[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = getPeak(vecRate,sRate.vecT,vecRestrictRange);
		sRate.dblPeakRate = dblPeakRate;
		sRate.dblPeakTime = dblPeakTime;
		sRate.dblPeakWidth = dblPeakWidth;
		sRate.vecPeakStartStop = vecPeakStartStop;
		sRate.intPeakLoc = intPeakLoc;
		sRate.vecPeakStartStopIdx = vecPeakStartStopIdx;
		sRate.dblOnset = [nan];
		
		intZetaIdxRate = min(max(1,intZETALoc-1),numel(vecRate));
		intZetaIdxInvRate = min(max(1,intPeakLocInvSign-1),numel(vecRate));
		
		if ~isnan(dblPeakTime)
			%assign array data
			
			%get onset
			[dblOnset,dblOnsetVal] = getOnset(vecRate,sRate.vecT,dblPeakTime,vecRestrictRange);
			sRate.dblOnset = dblOnset;
			vecLatencies = [dblMaxDTime dblMaxDTimeInvSign dblPeakTime dblOnset];
			vecLatencyVals = [vecRate(intZetaIdxRate) vecRate(intZetaIdxInvRate) vecRate(intPeakLoc) dblOnsetVal];
		end
	end
	
	%% build output structure
	sZETA = struct;
	sZETA.dblZETA = dblZETA;
	sZETA.dblD = dblD;
	sZETA.dblP = dblZetaP;
	sZETA.dblZetaT = dblMaxDTime;
	sZETA.intZetaIdx = intZETALoc;
	if boolStopSupplied
		sZETA.dblMeanZ = dblMeanZ;
		sZETA.dblMeanP = dblMeanP;
		sZETA.vecMu_Dur = vecMu_Dur;
		sZETA.vecMu_Pre = vecMu_Pre;
	end
	sZETA.vecSpikeT = vecSpikeT;
	sZETA.vecD = vecRealDiff;
	sZETA.vecRealFrac = vecRealFrac;
	sZETA.vecRealFracLinear = vecRealFracLinear;
	sZETA.cellRandT = cellRandT;
	sZETA.cellRandDiff = cellRandDiff;
	
	sZETA.dblD_InvSign = dblD_InvSign;
	sZETA.dblZetaT_InvSign = dblMaxDTimeInvSign;
	sZETA.intZetaIdx_InvSign = intPeakLocInvSign;
	sZETA.dblUseMaxDur = dblUseMaxDur;
	sZETA.vecLatencies = vecLatencies;
	sZETA.vecLatencyVals = vecLatencyVals;
	
	%latencies
	sLatencies = struct;
	sLatencies.Onset = vecLatencies(4);
	sLatencies.Peak = vecLatencies(3);
	sLatencies.ZETA = vecLatencies(1);
	sLatencies.ZETA_InvSign = vecLatencies(2);
	
	%% plot?
	if intPlot > 0
		sZETA.vecHandles = makeZetaPlot(vecSpikeTimes,matEventTimes,sZETA,sRate,intPlot);
	end
end

