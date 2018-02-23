'''''''''''''''''''''''''''''''''''''''''''''
Fast transients measured by Fermi and Swift show similar forms:
- Narrow pointed gaussian: Small_Gaussian
- Broad round gaussian: Simple_Gaussian
- Delta peak with exponential waste : Exponential

See plots of these Curves in the folder 'data/Template_Forms'

Templates = y-data of fast transient describes by above form
---> function elem Small_Gaussian, Simple_Gaussian, Exponential
---> function(Num_slices,noise_percentaige):
     - Num_slices = Number time slices = Number data points (int)
     - noise_percentage = range of noise in percentage (int/float)
'''''''''''''''''''''''''''''''''''''''''''''

import numpy as np


def Gauss(x, a, x0, sigma, b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b


@np.vectorize
def Step_exp(t, tp, Fp, b, alpha):
    if t < tp:
        return float(0)  # uniform data types
    if t >= tp:
        return float(Fp*(t/tp-b)**(-alpha))


### Fitting parameters averaged over 5-17 transients
Mean_params_Simple = np.array([0.626635975768735, 12.860830849499434, 4.6052849965482245, 0.1331112404132909])
Mean_params_Small = np.array([0.8957302067729599, 11.612091524250664, 0.9982224038570281, 0.11782671354588838])
Mean_params_Exp = np.array([2.86633333, 0.8355, -1.43666667,  1.41788333])


### Broad round gaussian
def broad_gaussian(Num_slices, noise_percentage):
    Timestamp = np.linspace(0, Num_slices, Num_slices)
    x_Gauss = range(int(4*Mean_params_Simple[1])) - Mean_params_Simple[1]
    # intepolation to guarantee correct form
    x_p = np.linspace(Timestamp[0], Num_slices, len(x_Gauss))
    y_p = (Gauss(x_Gauss, Mean_params_Simple[0], Mean_params_Simple[1], Mean_params_Simple[2], Mean_params_Simple[3]) - Mean_params_Simple[3]) / Mean_params_Simple[0] ## Normalized 1
    y = np.interp(Timestamp, x_p, y_p)
    # add noise
    noise = np.random.normal(0, noise_percentage/100.0, y.shape)
    y = y + noise
    return abs(y)


### small pointed gaussian
def narrow_gaussian(Num_slices, noise_percentage):
    Timestamp = np.linspace(0, Num_slices, Num_slices)
    x_Gauss = range(int(4*Mean_params_Small[1]))-Mean_params_Small[1]
    x_p = np.linspace(Timestamp[0], Num_slices, len(x_Gauss))
    y_p = (Gauss(x_Gauss, Mean_params_Small[0], Mean_params_Small[1], Mean_params_Small[2], Mean_params_Small[3])-Mean_params_Small[3])/Mean_params_Small[0]
    y = np.interp(Timestamp, x_p, y_p)
    noise = np.random.normal(0, noise_percentage/100.0, y.shape)
    y = y + noise
    return abs(y)


### exponential decay after delta Peak
def deltapeak_exponential(Num_slices, noise_percentage):
    Timestamp = np.linspace(0, Num_slices, Num_slices)
    x_Exp = np.linspace(-30, 170)
    x_p = np.linspace(Timestamp[0], Num_slices, len(x_Exp))
    y_p = Step_exp(x_Exp, Mean_params_Exp[0], Mean_params_Exp[1], Mean_params_Exp[2], Mean_params_Exp[3]) / Step_exp(x_Exp, Mean_params_Exp[0], Mean_params_Exp[1], Mean_params_Exp[2], Mean_params_Exp[3]).max()
    y = np.interp(Timestamp, x_p, y_p)
    noise = np.random.normal(0, noise_percentage/100.0, y.shape)
    y = y + noise
    return abs(y)


'''''''''''''''''''''''''''
Calling up the functions
E = Exponential(20,4)
S = Small_Gaussian(20,4)
SI = Simple_Gaussian(20,5)
'''''''''''''''''''''''''''
