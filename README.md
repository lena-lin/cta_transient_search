# CTA Transient Search
## Basis idea
Use Wavelet Denoising to find transients in CTA DL3 data. This data consists of 2d histograms containing number of gammalike events per ra/dec bin. Each histogram represents a short duration of observation (around 30 seconds). In this case a steady source (crab), a number of transients and the cosmic ray background is simulated. Wavelet Denoising is used to clean the
images/histograms from background-noise (see animation below). From the cleaned imgages a trigger criterion is extracted and monitored over time to find objects with high variability in CTA data.

## How to simulate ? 
Makefile settings: 
- irf path - path to CTA's IRF simulations 
- n_transients = Number simulated transients
- slices_per_part = Number timesteps for each simulation step: 1 transient = 1 simulation without transient, 1 simulation with transient, one simulation without transinet = 3*slices_per_part for 1 transient 
- transient_template_filename = random if exponential or gaussian, else 0 or 1 for gaussian an 2 for exponential 
- threshold = threshold for the trigger criterion in a.u.
- Ra, Dec of the transient if position is not randomly drag, see simulate_cube.py with parameter **-p**

## Branches 
### Different_Threshold_Studies
same setting compared to branch master with additional functions to plot some settings 
### new_background_studoes 
Additional component if a transient is simulated or not --> Easy way to get a false alarm rate 
### threshold studies 
Old branch for same studies as in Different_Threshold_Studies 
### Wavelet_studies
Change wavelet thresholding in k index and stationary methods 


## Additional Code
Simulations for finding usefull templates in different Repo stored [here](https://github.com/JanaMo/Transient_Simulation_Methods)

![wavelet](https://raw.githubusercontent.com/mackaiver/wavelet-denoising/sliding_bg_window/transient_sw.gif)
