# CTA Transient Search

Use Wavelet Denoising to find transients in CTA DL3 data. This data consists of 2d histograms containing number of gammalike events per ra/dec bin. Each histogram represents a short duration of observation (around 30 seconds). In this case a steady source (crab), a number of transients and the cosmic ray background is simulated. Wavelet Denoising is used to clean the
images/histograms from background-noise. From the cleaned imgages a trigger criterion is extracted and monitored over time to find objects with high variability in CTA data.

![wavelet](https://raw.githubusercontent.com/mackaiver/wavelet-denoising/sliding_bg_window/transient_sw.gif)
