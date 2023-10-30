function [vecTime,vecRate,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,boolUseParallel)
	%getIFR Returns instaneous firing rate. Syntax:
	%   [vecTime,vecRate,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,boolUseParallel)
	%Required input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: min of trial start to trial start]
	%
	%Optional inputs:
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of bins) [default: 2]
	%	- dblMinScale: minimum derivative scale in log-seconds [default: round(log(1/1000) / log(dblBase))]
	%	- dblBase: critical value for locally dynamic derivative [default: 1.5]
	%
	%Outputs:
	%	- vecTime; Time points corresponding to rates in vecRate
	%	- vecRate; Instantaneous firing rate in Hz
	%	- sIFR; structure with fields:
	%		- vecRate;
	%		- vecTime;
	%		- vecDiff;
	%		- vecScale; 
	%
	%v1.8 - 23 October 2023
	
	%Version history:
	%1.0 - 24 January 2019
	%	Created by Jorrit Montijn - split from getMultiScaleDeriv.m
	%1.1 - 24 June 2020
	%	Syntax cleanup [by JM]
	%1.2 - 3 July 2020
	%	Conform to ZETA naming [by JM]
	%1.3 - 14 October 2021
	%	Fixed bug [by JM]
	%1.4 - 16 February 2023
	%	Correction: updated function description to real default values [by JM]
	%1.5 - 26 May 2023
	%	Faster computation time, changed default parallel-processing behaviour [by JM]
	%1.6 - 6 June 2023
	%	Fixed mismatch of vecTime/vecRate after last update [by JM]
	%1.7 - 22 August 2023
	%	Changed default dblUseMaxDur to min instead of median ITI to match zetatest default  [by JM]
	%1.8 - 23 October 2023
	%	Removed plotting [by JM]
	
	%% set default values
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 2;
	end
	if ~exist('dblBase','var') || isempty(dblBase)
		dblBase = 1.5;
	end
	if ~exist('dblMinScale','var') || isempty(dblMinScale)
		dblMinScale = round(log(1/1000) / log(dblBase));
	end
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min(diff(vecEventStarts(:,1)));
	end
	if ~exist('boolUseParallel','var') || isempty(boolUseParallel)
		objPool = gcp('nocreate');
		if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
			boolUseParallel = false;
		else
			boolUseParallel = true;
		end
	end
	if size(vecEventStarts,2) > 2
		vecEventStarts = vecEventStarts';
	end
	
	%% get difference from uniform
	[vecRealDiff,vecThisSpikeFracs,vecThisFracLinear,vecTime] = ...
		getTempOffsetOne(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
	intSpikes = numel(vecTime);
	
	%% assign dummy output
	if numel(vecRealDiff) < 3
		vecTime = [];
		vecRate = [];
		sIFR = struct;
		sIFR.vecRate = vecRate;
		sIFR.vecTime = vecTime;
		sIFR.vecDiff = vecRealDiff;
		sIFR.vecScale = [];
		return
	end
	
	%% get multi-scale derivative
	intMaxRep = size(vecEventStarts,1);
	dblMeanRate = (intSpikes/(dblUseMaxDur*intMaxRep));
	[vecRate,sMSD] = getMultiScaleDeriv(vecTime,vecRealDiff,intSmoothSd,dblMinScale,dblBase,dblMeanRate,dblUseMaxDur,boolUseParallel);
	vecTime = sMSD.vecT;
	
	%% build output
	if nargout > 1
		sIFR = struct;
		sIFR.vecRate = vecRate;
		sIFR.vecTime = vecTime;
		sIFR.vecDiff = vecRealDiff;
		sIFR.vecScale = sMSD.vecScale;
	end
end

