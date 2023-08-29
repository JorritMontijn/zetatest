%runExampleZETA Run example ZETA-test
%
%This code loads data from an example cell and performs some analyses as a tutorial

%Version history:
%1.0 - 15 June 2020
%	Created by Jorrit Montijn
%1.1 - 23 Sept 2021
%	Added switch to use empirical null distribution for significance calculation [by JM]
%2.0 - 10 January 2022
%	Updated examples to conform to new syntax and added two-sample zeta [by JM]
%3.0 - 29 August 2023
%   Updated to match python example

%% note on performance
%{
The code has been computationally optimized, but since the ZETA-test is based on bootstraps, it
might still take a couple of seconds to compute. In cases where you recorded a lot of spikes, or a
long time-series, parallel processing can usually improve the speed. Unfortunately, this is not
*always* the case, so you might have to see whether parallel processing improves things for you on a
case-by-case basis. Note that due to extensive vectorization, MATLAB always outperforms Python.

In the examples below, I recorded the following processing times on my PC with default parameters:

MATLAB with parallel processing:
    zeta-test: 0.56 s
    time-series zeta-test: 2.77 s 

MATLAB without parallel processing:
    zeta-test: 1.09 s
    time-series zeta-test: 19.88 s

Python without parallel processing:
    zeta-test: 1.41 s
    time-series zeta-test: 32.24 s

%}

%% load data for example cell
rng(1,'twister'); %to make the output deterministic
sLoad = load('ExampleDataZetaTest.mat'); %loads matlab data file

%retrieve the spike times for as a vector from the field in sNeuron for two neurons
vecSpikeTimes1 = sLoad.sNeuron(1).SpikeTimes;
vecSpikeTimes2 = sLoad.sNeuron(2).SpikeTimes;

% load stimulation information
sStim = sLoad.sStim;
vecStimulusStartTimes = sStim.StimOnTime(:); %use (:) to ensure it's a column vector
vecStimulusStopTimes = sStim.StimOffTime(:);
matEventTimes = cat(2,vecStimulusStartTimes,vecStimulusStopTimes); % put stimulus start and stop times together into a [T x 2] matrix

%% calculate instantaneous firing rate without performing the ZETA-test
%if we simply want to plot the neuron's response, we can use:
[vecTime,vecRate,sIFR] = getIFR(vecSpikeTimes1,vecStimulusStartTimes);

%% run the ZETA-test with default parameters
%if we simply want to know if the neuron responds, no hassle, we can
%use this simple syntax with default parameters:
hTic=tic;
dblZetaP_default = zetatest(vecSpikeTimes1,vecStimulusStartTimes);
dblElapsedTime=toc(hTic);
fprintf('\nDefault parameters (elapsed time: %.2f s):\nzeta-test p-value: %f\n',dblElapsedTime,dblZetaP_default)

%% run the ZETA-test with specified parameters
%however, we can also specify the parameters ourselves
rng(1,'twister'); %to make the output deterministic
dblUseMaxDur = min(diff(vecStimulusStartTimes)); %minimum of trial-to-trial durations
intResampNum = 50; %~50 random resamplings should give us a good enough idea if this cell is responsive, but if it's close to 0.05, we should increase this #. Generally speaking, more is better, so let's put 100 here.
intPlot = 3;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
intLatencyPeaks = 4; %how many latencies do we want? 1=ZETA, 2=-ZETA, 3=peak, 4=first crossing of peak half-height
vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
boolDirectQuantile = false;%if true; uses the empirical null distribution rather than the Gumbel approximation. Note that in this case the accuracy of your p-value is limited by the # of resamplings
dblJitterSize = 1; %scalar value, sets the temporal jitter window relative to dblUseMaxDur (default: 2). Note that a value of "2" therefore means that the last data point that can be used is at t=last_onset + dblUseMaxDur*3

%then run ZETA with those parameters
hTic2=tic;
[dblZetaP,sZETA,sRate,vecLatencies] = zetatest(vecSpikeTimes1,matEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize);
dblElapsedTime2=toc(hTic2);
fprintf("\nSpecified parameters (elapsed time: %.2f s):\nzeta-test p-value: %f\nt-test p-value:%f\n",dblElapsedTime2,dblZetaP,sZETA.dblMeanP)

% Note on the latencies: while the peaks of ZETA and -ZETA can be useful for diagnostic purposes,
% they are difficult to interpret, so we suggest sticking to the peak time (vecLatencies[2]),
% which is more easily interpretable. Please read the paper for more information.

%% by popular demand: using a baseline that precedes the onset
%putting the baseline before the stimulus can be done by simply subtracting
%the baseline duration from the event times:
dblBaselineDuration = 0.5;
matEventTimesWithPrecedingBaseline = matEventTimes - dblBaselineDuration;
rng(1,'twister'); %to make the output deterministic

%then run ZETA with the new times
[dblZetaP_pb,sZETA_pb,sRate_pb,vecLatencies_pb] = zetatest(vecSpikeTimes1,matEventTimesWithPrecedingBaseline,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize);

%% however, the zeta function of course won't be able to tell the difference, so all timings are off by 500 ms.
%here we change the figure labels/titles (you can ignore this bit if you're not using the figure)
dblBaselineDurationMs = dblBaselineDuration*1000;
drawnow;hFig = gcf;
for intPlotNr=1:numel(hFig.Children)
	%adjust x-ticks
	if contains(hFig.Children(intPlotNr).XLabel.String,'Time ')
		set(hFig.Children(intPlotNr),'xticklabel',cellfun(@(x) num2str(str2double(x)-dblBaselineDuration),get(hFig.Children(intPlotNr),'xticklabel'),'UniformOutput',false));
	end
	%adjust timings in title
	strTitle = hFig.Children(intPlotNr).Title.String;
	[vecStart,vecStop]=regexp(strTitle,'[=].*?[m][s]');
	for intEntry=1:numel(vecStart)
		strOldNumber=hFig.Children(intPlotNr).Title.String((vecStart(intEntry)+1):(vecStop(intEntry)-2));
		strTitle = strrep(strTitle,strcat('=',strOldNumber,'ms'),strcat('=',num2str(str2double(strOldNumber)-dblBaselineDurationMs),'ms'));
	end
	hFig.Children(intPlotNr).Title.String = strTitle;
end
drawnow;

%here we adjust the times in the variables that getZeta returns
vecLatencies_pb = vecLatencies_pb - dblBaselineDuration;
sZETA_pb.vecSpikeT = sZETA_pb.vecSpikeT - dblBaselineDuration;
sRate_pb.vecT = sRate_pb.vecT - dblBaselineDuration;
sRate_pb.dblPeakTime = sRate_pb.dblPeakTime - dblBaselineDuration;
sRate_pb.dblOnset = sRate_pb.dblOnset - dblBaselineDuration;

%% run the time-series ZETA-test
% take subselection of data
intUseTrialNum = 480;
vecStimulusStartTimesTs = vecStimulusStartTimes(1:intUseTrialNum);
vecStimulusStopTimesTs = vecStimulusStopTimes(1:intUseTrialNum);
matEventTimesTs = cat(2,vecStimulusStartTimesTs,vecStimulusStopTimesTs);

% first transform the data to time-series
fprintf('\nRunning time-series zeta-test; This will take a couple of seconds\n')
dblStartT = 0;
dblEndT = vecStimulusStopTimesTs(end) + dblUseMaxDur*5;
dblSamplingRate = 50.0; % simulate acquisition rate
dblSampleDur = 1/dblSamplingRate;
vecTimestamps = dblStartT:dblSampleDur:dblEndT+dblSampleDur;
vecSpikesBinned = histcounts(vecSpikeTimes1, vecTimestamps);
vecTimestamps = vecTimestamps(1:(end-1));
dblSmoothSd = 1;
intSmoothRange = 2*ceil(dblSmoothSd);
vecFilt = normpdf(-intSmoothRange:intSmoothRange, 0, dblSmoothSd);
vecFilt = vecFilt / sum(vecFilt);

% pad array
intPadSize = floor(numel(vecFilt)/2);
vecData = padarray(vecSpikesBinned, [0 intPadSize],'replicate');

% filter
vecData = conv(vecData, vecFilt, 'valid');

%set random seed
rng(1,'twister');

% time-series zeta-test with default parameters
hTic3 = tic;
dblTsZetaP = zetatstest(vecTimestamps, vecData, vecStimulusStartTimesTs);
dblElapsedTime3=toc(hTic3);
fprintf("\nDefault parameters (elapsed time: %.2f s):\ntime-series zeta-test p-value: %f\n",dblElapsedTime3,dblTsZetaP)

%% run time-series zeta-test with specified parameters
% set random seed
rng(1,'twister');
fprintf('\nRunning time-series zeta-test with specified parameters; This will take a couple of seconds\n')

% run test
hTic4 = tic;
[dblTsZetaP2, sZetaTs] = zetatstest(vecTimestamps, vecData, matEventTimesTs,dblUseMaxDur, intResampNum, intPlot, boolDirectQuantile, dblJitterSize);
dblElapsedTime4 = toc(hTic4);
fprintf("\nSpecified parameters (elapsed time: %.2f s):\ntime-series zeta-test p-value: %f\nt-test p-value: %f\n",dblElapsedTime4,dblTsZetaP,sZetaTs.dblMeanP)

return
%% finally, run the two-sample ZETA-test
% NOTE: THIS IS STILL UNDER DEVELOPMENT! USE AT YOUR OWN RISK
%case 1: are neurons 1 & 2 responding differently to a set of visual stimuli?
%in this example, we are testing two simultaneously recorded neurons, so we can use the paired
%version of the test (see note after typing "help zetatest2"). However, since these neurons are not
%very noisy, paired testing should have little effect on the statistical significance.
intResampNum = 500; %the two-sample test is more variable, as it depends on differences, so it requires more resamplings
boolPairedTest = true;
[dblZetaPaired,sZETA2] = zetatest2(vecSpikeTimes1,matEventTimes,vecSpikeTimes2,matEventTimes,boolPairedTest,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize);

%case 2a: is neuron 1 responding differently to gratings oriented at 30 and 60 degrees?
intPlot = 3;
boolPairedTest = false;
vecTrials1 = sStim.Orientation==30;
vecTrials2 = sStim.Orientation==60;
dblZetaUnpaired1 = zetatest2(vecSpikeTimes1,matEventTimes(vecTrials1,:),vecSpikeTimes1,matEventTimes(vecTrials2,:),boolPairedTest,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize);

%case 2b: is neuron 2 responding differently to gratings oriented at 30 and 60 degrees?
dblZetaUnpaired2 = zetatest2(vecSpikeTimes2,matEventTimes(vecTrials1,:),vecSpikeTimes2,matEventTimes(vecTrials2,:),boolPairedTest,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize);
