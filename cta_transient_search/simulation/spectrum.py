import numpy as np
import astropy.units as u
from scipy.integrate import quad
from ebltable.tau_from_model import OptDepth

#tau = OptDepth.readmodel(model = 'dominguez')


def random_power(e_min, e_max, N):
    '''
    Returns random numbers from power law with given index and energy range (e_min, e_max in TeV)
    '''
    index = 2.48
    u = np.random.uniform(0, 1, N)
    # inverse transform sampling
    return ((e_max**(1-index) - e_min**(1-index)) * u + e_min**(1-index))**(1/(1-index))


def number_particles_crab(T, e_min, e_max, simulation_area):
    '''
    Returns the number of particles arriving from a pointlike source with known power law distribution

    Parameters:
    T: observation time
    e_min: lower energy bound
    e_max: upper energy bound
    simulation_area: simulated area (twice CTA radius)
    C: Constant flux factor
    index: power lax index
    '''

    # Spectrum HEGRA (2004)
    # index = 2.62
    # C = 2.83e-11 / (u.cm**2 * u.s * u.TeV)
    # E_0 = 1 * u.TeV

    # Spectrum used in ctools, energy in MeV!!!!!!!!
    index = 2.48
    C = 5.7e-16 / (u.cm**2 * u.s * u.MeV)
    E_0 = 0.3e6 * u.MeV
    e_min = e_min.to(u.MeV)
    e_max = e_max.to(u.MeV)

    return int((simulation_area.to(u.cm**2)*C*T*E_0**(index)/(1-index))*((e_max)**(1-index) - (e_min)**(1-index)))


def number_particles_trans(T, e_min, e_max, simulation_area, Redshift):
    '''
    Returns the number of particles arriving from a pointlike source with known power law distribution and EBL attenuation
    Parameters:
    T: observation time
    e_min: lower energy bound
    e_max: upper energy bound
    simulation_area: simulated area (twice CTA radius)
    C: Constant flux factor
    index: power lax index
    z: Redshift
    '''
    index = 2.48
    C = 5.7e-16 / (u.cm**2 * u.s * u.MeV)
    E_0 = 0.3e6 * u.MeV
    e_min = e_min.to(u.MeV)
    e_max = e_max.to(u.MeV)

    simulated_area = simulation_area.to(u.cm**2)
    sim_area = simulated_area.value
    mean_e = (emax+emin)/2
    z = Redshift[0]

    '''
    integrate with quad to get more differenced EBL attenuation:
    integrate:  E**(-index)*np.exp(-1.0*tau.opt_depth(z,E)) from E_min to E_max
    '''
    def spectrum_crab_ebl(z):
         return lambda E : (C * (E/E_0) ** (-index) * np.exp(-1.0*tau.opt_depth(z,E.value))).value

    I = quad(spectrum_crab_ebl(z),emin,emax)[0]

    print('Simulated EBL attenuation with z= ',z,
    '\n N without EBL =',int((sim_area*C*T.value*E_0**(index)/(1-index))*((emax)**(1-index) - (emin)**(1-index))),
    '\n N is now = ',int((sim_area*C*T.value*E_0**(index))*I))
    return int((sim_area*C*T.value*E_0**(index))*I)
