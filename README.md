# The ZETA-test repository
Repository containing ZETA-test functions and dependencies. For an example of how to use the code, check runExampleZETA.m in the /example/ subfolder. Your output should look like the images in the same directory. 

The article describing the original ZETA-test has been published in eLife: https://elifesciences.org/articles/71969

If you're looking for the original ZETA repository, you can find it here: https://github.com/JorritMontijn/ZETA

Note on updates and maintenance: the functions provided in this repository are in different stages of development: the zetatest has been extensively tested on real and artificial data, and has been peer-reviewed, so we are confident whatever it produces is reliable and statistically sound. Compared to the getZeta.m in the original repository, this version has a somewhat different syntax, and adds a data-stitching step to ensure that time is contiguous. If your data was already contiguous, you should notice no difference, but if your events are more randomly distributed, you might get better performance with this version. 
The time-series zetatstest has also been tested on a variety of artificial and real benchmarks and seems to perform better than t-tests and one-way ANOVAs, so we are also fairly confident of the zetatstest, but this method has not yet been peer-reviewed (preprint will be published before the summer). The two-sample tests are still under active development, so you are welcome to try them out and let us know what you find of their performance, but we make no promises as to their performance - it might well be possible we'll scrap them altogether, so use these two-sample tests at your own risk. More information on these tests can be found in runExampleZETA.m and the help comments of the respective functions.

 
This repository contains five main functions:
1) zetatest.m: Calculates the Zenith of Event-based Time-locked Anomalies (ZETA) for spike times of a single neuron. Outputs a p-value.
2) zetatstest.m: Calculates the time-series version of ZETA, for data such as calcium imaging or EEG recordings.
3) zetatest2.m: Same as (1), but for testing whether two neurons respond differently to the same stimulus; or whether one neuron responds differently to two sets of stimuli. Still under construction.
4) zetatstest2.m: Same as (3), but for testing differences between two time-series data arrays. Still under construction.
5) getIFR.m: Calculates the instantaneous firing rate (IFR) without running the ZETA-test. Use this as you would a PSTH function.

# Rationale for ZETA

Neurophysiological studies depend on a reliable quantification of whether and when a neuron responds to stimulation, be it sensory, optogenetically or otherwise. However, current statistical analysis methods to determine a neuron’s responsiveness require arbitrary parameter choices, such as a binning size. This choice can change the results of the analysis, which invites bad statistical practice and reduces the replicability of analyses. Moreover, many methods, such as bin-wise t-tests, only detect classically mean-rate modulated  cells. Especially with advent of techniques that yield increasingly large numbers of cells, such as Neuropixels  recordings , it is important to use tests for cell-inclusion that require no manual curation. Here, we present the parameter-free ZETA-test, which outperforms common approaches, in the sense that it includes more cells at a similar false-positive rate. 
Finally, ZETA’s timescale-, parameter- and binning-free nature allowed us to implement a ZETA-derived algorithm (using multi-scale derivatives) to calculate peak onset and offset latencies in neuronal spike trains with theoretically unlimited temporal resolution. 

We have now adapted this approach to time-series data as well, and our initial benchmarks are very promising. That said, the zetatstest and zetatstest2 are still undergoing testing, so we are not yet making any claims as to their performance or correctness. Use at your own risk.

Please send any questions or comments to j.montijn at nin.knaw.nl.


# Dependencies
The ZETA-test functions require the following Mathworks toolboxes to work:
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- (Optional: Parallel Computing Toolbox to reduce computation time)


![zeta_image](https://user-images.githubusercontent.com/15422591/135059690-2d7f216a-726e-4080-a4ec-2b3fae78e10c.png)
