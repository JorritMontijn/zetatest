function [dblOnset,dblValue,dblBaseVal,dblPeakT,dblPeakVal] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange)
	%getOnset Returns peak onset. Syntax:
	%    [dblOnset,dblValue,dblBaseVal,dblPeakT,dblPeakVal] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange)
	%
	%Required input:
	%	- vecData [N x 1]: values
	%
	%Optional inputs:
	%	- vecT [N x 1]: timestamps corresponding to vecData (default: [1:N])
	%	- dblPeakT (float): timestamp corresponding to peak
	%	- vecRestrictRange [2 x 1]: restrict peak to lie within vecRestrictRange(1) and vecRestrictRange(end)
	%
	%Outputs:
	%	- dblOnset: time of peak onset (first crossing half-height of peak)
	%	- dblValue: value at peak onset
	%
	%Version history:
	%1.0 - February 26 2020
	%	Created by Jorrit Montijn
	%2.0 - August 22 2023
    %   Now actually uses dblPeakT [by JM]

	%%
	%check inputs
	if ~exist('vecT','var') || isempty(vecT)
		vecT = 1:numel(vecData);
	end
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [min(vecT) min(vecT)+range(vecT)];
	end
	
	%remove time points outside restricted range
	indRemove = vecT < vecRestrictRange(1) | vecT > vecRestrictRange(end);
	vecCropT = vecT(~indRemove);
	vecDataCropped = vecData(~indRemove);
	
	%find peak if none supplied
    if ~exist('dblPeakT','var') || isempty(dblPeakT) 
        [dblPeakVal,intPeakT] = max(vecDataCropped);
	    dblPeakT = vecT(intPeakIdx);
    else
        intPeakT = find(vecCropT > dblPeakT,1,'first');
    end
    dblPeakVal = vecDataCropped(intPeakT);
	
    %calculate first timepoint crossing half-height of peak 
	dblBaseVal = vecDataCropped(1);
	dblThresh = (dblPeakVal - dblBaseVal)/2 + dblBaseVal;
	if dblThresh > 0
		intOnsetIdx = find(vecDataCropped >= dblThresh,1,'first');
	else
		intOnsetIdx = find(vecDataCropped <= dblThresh,1,'first');
	end
	dblOnset = vecCropT(intOnsetIdx);
	dblValue = vecDataCropped(intOnsetIdx);
	
	%check if empty
	if isempty(intOnsetIdx)
		dblOnset = nan;
		dblValue = nan;
	end
end
