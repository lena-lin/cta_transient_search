# CTA Transient Search - New Templates 

New working branch for integration of new templates (light curves and spectra) and extension of time scales :
- Installation of functions for exponential and exponential light curves (See Templates in */Plots*)
- Variable setting of the transient start time depending on the selected template 
- Randomized selection of the template / the light curve shape 
- Check influence of noise 
- Comparison of Crab and transient peak fluxes --> Restriction of Crab Unit Flare(See examples 2 CU and 0.3 CU in */Plots*)
- Lowering and Changing of alter threshold 
- Tests : distance transient - steady source, time distance of tweo transients, strength of the transient in CUs, Threshold, sliding window gap.... 


<br><br>
Master branch from lena-lin: <br>
Use Wavelet Denoising to find transients in CTA DL3 data. This data consists of 2d histograms containing number of gammalike events per ra/dec bin. Each histogram represents a short duration of observation (around 30 seconds). In this case a steady source (crab), a number of transients and the cosmic ray background is simulated. Wavelet Denoising is used to clean the
images/histograms from background-noise. From the cleaned imgages a trigger criterion is extracted and monitored over time to find objects with high variability in CTA data.

![wavelet](https://raw.githubusercontent.com/mackaiver/wavelet-denoising/sliding_bg_window/transient_sw.gif)
