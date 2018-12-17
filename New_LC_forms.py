import numpy as np
from random import uniform

'''''''''''''''''''''''''''''''''''''''''''''
Fast transients measured by Fermi and Swift show similar forms:
- Narrow pointed gaussian: Small_Gaussian
- Broad round gaussian: Simple_Gaussian
- Delta peak with exponential waste : Exponential

See plots of these Curves in the folder 'data/Template_Forms'

Templates = y-data of fast transient describes by above form

---> function simulate_Exponential / simulate_Gaussian
---> simulate Gaussian(sigma_range_min,sigma_range_max,Num_slices,time_per_slice):
     - Num_slices = Number time slices = Number data points (int)
     - time_per_slice = time in seconds per slice
     - sigma_range_min = lowest sigma value taken from fits to Fermi data
     - sigma_range_max = highest sigma value taken from fits to Fermi data
---> simulate_Exponential(alpha_min,alpha_max,b_min,b_max,Num_slices,time_per_slice):
     - alpha_min = Minimum alpha value for ideal exponential decay within time window
     - alpha_max = Maximum alpha value for ideal exponential decay within time window
     - b_min = Minimum Shift in exponential decay taken from fits to Fermi data
     - b_max = Maximum Shift in exponential decay taken from fits to Fermi data
'''''''''''''''''''''''''''''''''''''''''''''


def Gauss_Norm(x, x0, sigma, b):
    return np.exp(-(x-x0)**2/(2*sigma**2))+b


def x_from_y_Gauss(percentage, y, sigma, mu):
    return -np.sqrt(2*np.log(1/(percentage*y)))*sigma+mu


@np.vectorize
def Step_exp(t, tp, Fp, b, alpha):
    if t < tp:
        return float(0)  # uniform data types
    if t >= tp:
        return float(Fp*(t/tp+b)**(-alpha))


def simulate_Gaussians(sigma_range_min, sigma_range_max, Num_slices, time_per_slice):
    '''
    Prefered start values:
    1.8 < sigma < 16 for brighter simulate_Gaussians
    0.45 < sigma < 2.2 for small Gaussians
    Num_slices = 40
    time_per_slice = 10
    '''
    time = Num_slices*time_per_slice

    '''
    Rescale sigmas to given time & randomply drag a float sigma value
    '''
    sigma_range_min_time = sigma_range_min/(600)*time
    sigma_range_max_time = sigma_range_max/(600)*time

    sigma = uniform(sigma_range_min_time, sigma_range_max_time)

    '''
    Adjust possible positions for the peak: mu
    Start with 3 slices + 1 sigma before peak
    End with 3 sliced and 1 sigma after peak
    '''
    mu_min = sigma_range_max_time+time_per_slice*2
    mu_max = Num_slices * time_per_slice - 2 * time_per_slice - sigma_range_max_time

    mu = uniform(mu_min, mu_max)

    '''
    Calculate #N_slices values, normal distributed around random mu and sigma
    '''
    xlin = np.linspace(0, time, Num_slices)
    ylin = Gauss_Norm(xlin, mu, sigma, 0)
    y = ylin.max()
    ylin = ylin/y
    '''
    Calculate True Startslice of the transient
    '''
    True_Start = int(x_from_y_Gauss(0.8, y, sigma, mu)/time_per_slice)  # Number of slice
    return ylin, True_Start


def simulate_Exponential(alpha_min, alpha_max, b_min, b_max, Num_slices, time_per_slice):
    '''
    Prefered start values:
    3 < alpha < 6
    0 < b < 2
    '''
    time = Num_slices*time_per_slice
    alpha_min_time = alpha_min
    b_min_time = b_min/600*time
    alpha_max_time = alpha_max
    b_max_time = b_max/600*time

    b = uniform(b_min_time, b_max_time)
    alpha = uniform(alpha_min_time, alpha_max_time)

    tp = uniform(0, time/3)
    params = np.array([tp, 1, b, alpha])

    x2 = np.linspace(0, time, Num_slices)
    y2 = Step_exp(x2, params[0], params[1], params[2], params[3])
    y2 = y2/y2.max()
    True_Start = int(tp/time_per_slice)+1

    return y2, True_Start
