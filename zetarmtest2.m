function [dblZetaP,sZETA] = zetarmtest2(vecT1,matCond1,vecT2,matCond2,...
		intResampNum,intPlot,boolDirectQuantile)
	%zetatstest2 Calculates difference in responsiveness index zeta for two timeseries
	%syntax: [dblZetaP,sZETA] = zetarmtest2(varT1,matCond1,varT2,matCond2,...
	%	intResampNum,intPlot,boolDirectQuantile)	

	%% prep data
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
	
	%% check input
	assert(ndims(matCond1) == 2 && ndims(matCond2) == 2 ...
		&& size(vecT1,2) == size(matCond1,2) && size(vecT2,2) == size(matCond2,2),...
		[mfilename ':InputError'],'Input sizes do not match');
	
	%% get ts-zeta diff
	[vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZetaIdx,matTracePerTrial1,matTracePerTrial2] = ...
		calcRmZetaTwo(vecT1,matCond1,vecT2,matCond2,intResampNum,boolDirectQuantile);
	
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
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = 0;
		sZETA.dblP = 1;
		sZETA.dblPeakT = nan;
		sZETA.intPeakIdx = [];
		if boolStopSupplied
			sZETA.dblMeanD = 0;
			sZETA.dblMeanP = 1;
		end
		sZETA.vecRefT = [];
		sZETA.vecD = [];
		sZETA.matRandD = [];
		
		sZETA.dblD_InvSign = 0;
		sZETA.dblPeakT_InvSign = nan;
		sZETA.intPeakIdx_InvSign = [];
		sZETA.dblUseMaxDur = nan;
		return
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
			imagesc(vecRefT,1:size(matTracePerTrial1,1),matTracePerTrial1,[0 1]);
			colormap(hot);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Rescaled data, condition 1');
			colorbar;
			grid off;
			
			subplot(2,3,4)
			imagesc(vecRefT,1:size(matTracePerTrial2,1),matTracePerTrial2,[0 1]);
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
		xlim([min(vecRefT) max(vecRefT)]);
		
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
		title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		xlim([min(vecRefT) max(vecRefT)]);
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
		sZETA.vecRefT = vecRefT;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
	end
end
