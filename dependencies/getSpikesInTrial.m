function [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts,dblMaxDur,boolReturnCells)
	%getSpikesInTrial Retrieves spiking times per trial
	%syntax: [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts,dblMaxDur,boolReturnCells)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%	- dblTrialDur: (optional) if supplied, uses trial duration for
	%			computations instead of assigning spikes directly to trials
	%	- boolReturnCells: (optional, default=false) if true, returns cell arrays for a single neuron 
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	%2.0 - February 7 2020
	%	Added overlap by using dblMaxDur [by JM]
	%2.1 - 24 May 2023
	%	Added support for cell array of spikes as input [by JM]
	
	%check inputs
	if nargin < 3
		dblMaxDur = [];
	end
	if ~exist('boolReturnCells','var') || isempty(boolReturnCells)
		boolReturnCells = false;
	end
	
	if iscell(vecSpikes)
		%run recurrently
		vecTrialPerSpike = cell(size(vecSpikes));
		vecTimePerSpike = cell(size(vecSpikes));
		for intCell=1:numel(vecSpikes)
			[vecTrialPerSpike{intCell},vecTimePerSpike{intCell}] = getSpikesInTrial(vecSpikes{intCell},vecTrialStarts,dblMaxDur);
		end
	else
		if ~isempty(dblMaxDur)
			%sort spikes
			intTrials = numel(vecTrialStarts);
			cellTrialPerSpike = cell(size(vecTrialStarts));
			cellTimePerSpike = cell(size(vecTrialStarts));
			for intTrial=1:intTrials
				%get spikes
				vecTheseSpikes = vecSpikes((vecSpikes >= vecTrialStarts(intTrial)) & vecSpikes < (vecTrialStarts(intTrial)+ dblMaxDur));
				vecTheseSpikes = vecTheseSpikes - vecTrialStarts(intTrial);
				
				% assign
				cellTrialPerSpike{intTrial} = intTrial*ones(size(vecTheseSpikes));
				cellTimePerSpike{intTrial} = vecTheseSpikes;
			end
			%transform to vectors & reorder
			if boolReturnCells
				vecTimePerSpike = cellTimePerSpike;
				vecTrialPerSpike = cellTrialPerSpike;
			else
				vecTimePerSpike = cell2vec(cellTimePerSpike);
				vecTrialPerSpike = cell2vec(cellTrialPerSpike);
			end
		else
			%sort spikes
			vecTrialPerSpike = nan(size(vecSpikes));
			vecTimePerSpike = nan(size(vecSpikes));
			for intSpike=1:numel(vecSpikes)
				%% build trial assignment
				vecTrialPerSpike(intSpike) = sum(vecTrialStarts < vecSpikes(intSpike));
				if vecTrialPerSpike(intSpike) > 0
					dblRemTime = vecTrialStarts(vecTrialPerSpike(intSpike));
				else
					dblRemTime = 0;
				end
				vecTimePerSpike(intSpike) = vecSpikes(intSpike) - dblRemTime;
			end
		end
	end
end

