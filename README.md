# ZETA-test
Repository containing ZETA-test functions and dependencies. For an example of how to use the code, check runExampleZETA.m in the /example/ subfolder. Your output should look like the images in the same directory. 

The article describing the original ZETA-test has been published in eLife: https://elifesciences.org/articles/71969

If you're looking for the original ZETA repository, you can find it here: https://github.com/JorritMontijn/ZETA

Note on updates and maintenance: we are currently still working on a time-series version of ZETA for calcium imaging data, so you might notice we are regularly updating the repository. zetatstest.m is still undergoing testing and regular changes, so we make no claims regarding its performance. zetatest.m, however, is stable and well tested, so you can safely use it. Note that the syntax is changed compared to the original getZeta.m file to make it conform to the other functions. In addition, this new repository has added two-sample versions of both the ZETA-test (zetatest2.m) and and the time-series ZETA-test (zetatstest2.m). More information on these tests can be found in runExampleZETA.m and the help comments of the respective functions.

 
This repository contains five main functions:
1) zetatest.m: Calculates the Zenith of Event-based Time-locked Anomalies (ZETA) for spike times of a single neuron. Outputs a p-value.
2) zetatest2.m: Same as above, but for testing whether two neurons respond differently to the same stimulus; or whether one neuron responds differently to two sets of stimuli.
3) zetatstest.m: Calculates the time-series version of ZETA, for data such as calcium imaging or EEG recordings.
4) zetatstest2.m: Same as above, but for testing differences between two time-series data arrays.
5) getIFR.m: Calculates the instantaneous firing rate (IFR) without running the ZETA-test. Use this as you would a PSTH function.

Rationale for ZETA

Neurophysiological studies depend on a reliable quantification of whether and when a neuron responds to stimulation, be it sensory, optogenetically or otherwise. However, current statistical analysis methods to determine a neuron’s responsiveness require arbitrary parameter choices, such as a binning size. This choice can change the results of the analysis, which invites bad statistical practice and reduces the replicability of analyses. Moreover, many methods, such as bin-wise t-tests, only detect classically mean-rate modulated  cells. Especially with advent of techniques that yield increasingly large numbers of cells, such as Neuropixels  recordings , it is important to use tests for cell-inclusion that require no manual curation. Here, we present the parameter-free ZETA-test, which outperforms common approaches, in the sense that it includes more cells at a similar false-positive rate. 
Finally, ZETA’s timescale-, parameter- and binning-free nature allowed us to implement a ZETA-derived algorithm (using multi-scale derivatives) to calculate peak onset and offset latencies in neuronal spike trains with theoretically unlimited temporal resolution. 

We have now adapted this approach to time-series data as well, and our initial benchmarks are very promising. That said, the zetatstest and zetatstest2 are still undergoing testing, so we are not yet making any claims as to their performance or correctness. Use at your own risk.

Please send any questions or comments to j.montijn at nin.knaw.nl.


Dependencies
The ZETA-test functions require the following Mathworks toolboxes to work:
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- (Optional: Parallel Computing Toolbox to reduce computation time)


![zeta_image](https://user-images.githubusercontent.com/15422591/135059690-2d7f216a-726e-4080-a4ec-2b3fae78e10c.png)
