# The ZETA-test repository
Repository containing ZETA-test functions and dependencies. Most up to date version is located at https://github.com/Herseninstituut/zetatest 

For an example of how to use the code, check runExampleZETA.m in the /example/ subfolder. Your output should look like the images in the same directory. 

A pre-print describing data-stitching, the time-series ZETA-test, and the two-sample tests is online: https://www.biorxiv.org/content/10.1101/2023.10.30.564780v1

The article describing the original ZETA-test has been published in eLife: https://elifesciences.org/articles/71969

If you're looking for the original ZETA repository, you can find it here: https://github.com/JorritMontijn/ZETA

The ZETA-test for spiking data has been extensively tested on real and artificial data, and has been peer-reviewed. The time-series ZETA-test and two-sample ZETA-tests are also thoroughly tested and are described in our pre-print, which we will submit for peer-review soon. We are confident all ZETA-tests, including the two-sample tests, are reliable and statistically sound, but if you find any bugs, do let us know on the Issues page!

More information on these tests can be found in runExampleZETA.m and the help comments of the respective functions.

This repository contains five main functions:
1) zetatest.m: Calculates the Zenith of Event-based Time-locked Anomalies (ZETA) for spike times of a single neuron. Outputs a p-value.
2) zetatstest.m: Calculates the time-series version of ZETA, for data such as calcium imaging or EEG recordings.
3) zetatest2.m: Same as (1), but for testing whether two neurons respond differently to the same stimulus; or whether one neuron responds differently to two sets of stimuli.
4) zetatstest2.m: Same as (2), but for testing differences between two time-series data arrays.
5) getIFR.m: Calculates the instantaneous firing rate (IFR) without running the ZETA-test. Use this as you would a PSTH function.

Additionally, we provide two ZETA-based functions for latency estimation:
1) latenzy.m: Estimates the response latency for spike times of a single neuron.
2) latenzy2.m: Estimates when spiking starts to diverge between two conditions.


# Rationale for ZETA

Neurophysiological studies depend on a reliable quantification of whether and when a neuron responds to stimulation, be it sensory, optogenetically or otherwise. However, current statistical analysis methods to determine a neuron’s responsiveness require arbitrary parameter choices, such as a binning size. This choice can change the results of the analysis, which invites bad statistical practice and reduces the replicability of analyses. Moreover, many methods, such as bin-wise t-tests, only detect classically mean-rate modulated  cells. Especially with advent of techniques that yield increasingly large numbers of cells, such as Neuropixels  recordings or two-photon calcium imaging, it is important to use tests for cell-inclusion that require no manual curation. Here, we present a new family of statistical tests for responses in point-event and time-series data for one- and two-sample comparisons: the family of ZETA-tests. As shown in our papers, they outperform approaches such as optimally-binned ANOVAs, t-tests and model-based approaches, in the sense that it includes more cells in real neurophysiological data at a similar false-positive rate. 

## Latency

For computing the latency of a response, we recommend to use the LatenZy test (https://github.com/Herseninstituut/latenZy), based on the zeta-test. 

## Credits

Original code by Jorrit Montijn. Maintenance by Alexander Heimel. Please cite https://elifesciences.org/articles/71969 if you use the test in a publication.


# Dependencies
The ZETA-test functions require the following Mathworks toolboxes to work:
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- (Optional: Parallel Computing Toolbox to reduce computation time)


![zeta_image](https://user-images.githubusercontent.com/15422591/135059690-2d7f216a-726e-4080-a4ec-2b3fae78e10c.png)
