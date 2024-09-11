function vecHandles = makeZetaPlot(vecSpikeTimes,matEventTimes,sZETA,sRate,intPlot)
	%makeZetaPlot Makes zeta plot
	%   vecHandles = makeZetaPlot(vecSpikeTimes,matEventTimes,sZETA,sRate,intPlot)
	
	%% check input
	if intPlot == 0,return;end
	vecHandles = nan(1,6);
	
	%% unpack
	dblZETA = sZETA.dblZETA;
	dblZetaD = sZETA.dblD;
	dblZetaP = sZETA.dblP;
	dblZetaT = sZETA.dblZetaT;
	intZetaIdx = sZETA.intZetaIdx;
	if isfield(sZETA,'dblMeanZ')
		boolStopSupplied = true;
		dblMeanZ = sZETA.dblMeanZ;
		dblMeanP = sZETA.dblMeanP;
		vecMu_Dur = sZETA.vecMu_Dur;
		vecMu_Pre = sZETA.vecMu_Pre;
	else
		boolStopSupplied = false;
	end
	vecSpikeT = sZETA.vecSpikeT;
	vecRealFrac = sZETA.vecRealFrac;
	vecRealFracLinear = sZETA.vecRealFracLinear;
	vecRealDiff = sZETA.vecD;
	cellRandT = sZETA.cellRandT;
	cellRandDiff = sZETA.cellRandDiff;
	
	dblZetaD_InvSign = sZETA.dblD_InvSign;
	dblZetaT_InvSign = sZETA.dblZetaT_InvSign;
	intZetaIdx_InvSign = sZETA.intZetaIdx_InvSign;
	dblUseMaxDur = sZETA.dblUseMaxDur;
	vecLatencies = sZETA.vecLatencies;
	vecLatencyVals = sZETA.vecLatencyVals;
	
	%% MSD
	if ~isempty(sRate)
		vecRate = sRate.vecRate;
		vecT = sRate.vecT;
		vecM = sRate.vecM;
		vecScale = sRate.vecScale;
		matMSD = sRate.matMSD;
		vecV = sRate.vecV;
	end
	
	%% plot
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
		hRaster = subplot(2,3,1);
		vecHandles(1) = hRaster;
		plotRaster(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,10000);
		xlabel('Time after event (s)');
		ylabel('Trial #');
		title('Spike raster plot');
		grid off;
	end
	
	%plot
	vecHandles(2) = subplot(2,3,2);
	sOpt = struct;
	sOpt.handleFig =-1;
	if dblUseMaxDur < 0.5
		dblBinSize = dblUseMaxDur/40;
	else
		dblBinSize = 0.025;
	end
	vecBins = 0:dblBinSize:dblUseMaxDur;
	warning('off','doPEP:WrongSyntax');
	[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,vecBins,matEventTimes(:,1),sOpt);
	warning('on','doPEP:WrongSyntax');
	errorbar(vecWindowBinCenters,vecMean,vecSEM);
	ylim([0 max(get(gca,'ylim'))]);
	title(sprintf('Mean spiking over trials'));
	xlabel('Time after event (s)');
	ylabel('Mean spiking rate (Hz)');
	
	vecHandles(3) = subplot(2,3,3);
	plot(vecSpikeT,vecRealFrac)
	hold on
	plot(vecSpikeT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
	title(sprintf('Real data'));
	xlabel('Time after event (s)');
	ylabel('Fractional position of spike in trial');
	
	vecHandles(4) = subplot(2,3,4);
	cla;
	hold all
	for intIter=1:intPlotIters
		plot(cellRandT{intIter},cellRandDiff{intIter},'Color',[0.5 0.5 0.5]);
	end
	plot(vecSpikeT,vecRealDiff,'Color',lines(1));
	scatter(dblZetaT,vecRealDiff(intZetaIdx),'bx');
	scatter(dblZetaT_InvSign,vecRealDiff(intZetaIdx_InvSign),'b*');
	hold off
	xlabel('Time after event (s)');
	ylabel('Offset of data from linear (s)');
	if boolStopSupplied
		title(sprintf('ZETA=%.3f (p=%.3f), z(Hz)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanZ,dblMeanP));
	else
		title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
	end
	ylim([-1.001 1.001]*max(abs(get(gca,'ylim'))));
	
	%% plot
	if intPlot == 1
		vecHandles(1) = subplot(2,3,1);
		stairs(vecT,vecRate)
		xlabel('Time after event (s)');
		ylabel('Instantaneous firing rate (Hz)');
		title(sprintf('Peri Event Plot (PEP)'));
	elseif intPlot > 1
		vecHandles(5) = subplot(2,3,5);
		imagesc(matMSD');
		set(gca,'ytick',[]);
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Timestamp index (#)');
		title('Multi-scale derivatives');
		grid off
		
		vecHandles(6) = subplot(2,3,6);
		if numel(vecT) > 10000
			vecSubset = round(linspace(1,numel(vecT),10000));
			plot(vecT(vecSubset),vecRate(vecSubset));
		else
			plot(vecT,vecRate);
		end
		xlabel('Time after event (s)');
		ylabel('Instantaneous firing rate (Hz)');
		title(sprintf('Peri Event Plot (PEP)'));
	end
	xlim([0 dblUseMaxDur]);
	
	hold on;
	scatter(vecLatencies(1),vecLatencyVals(1),'bx'); %ZETA
 	scatter(vecLatencies(2),vecLatencyVals(2),'b*'); %-ZETA
	scatter(vecLatencies(3),vecLatencyVals(3),'gx'); %Peak
	scatter(vecLatencies(4),vecLatencyVals(4),'rx'); %Onset
	title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,On=%.2fms',vecLatencies(1)*1000,vecLatencies(2)*1000,vecLatencies(3)*1000,vecLatencies(4)*1000));
	hold off
	xlim([0 dblUseMaxDur]);
	
	if intPlot > 3
		axes(hRaster);
		vecY = get(gca,'ylim');
		hold on;
		plot(vecLatencies(1)*[1 1],vecY,'b--');
		plot(vecLatencies(2)*[1 1],vecY,'b-.');
		plot(vecLatencies(3)*[1 1],vecY,'g--');
		plot(vecLatencies(4)*[1 1],vecY,'r--');
		hold off
	end
	
	fixfig;
end


