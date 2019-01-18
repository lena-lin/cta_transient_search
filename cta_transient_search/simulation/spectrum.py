import numpy as np
import astropy.units as u
from scipy.integrate import trapz
from ebltable.tau_from_model import OptDepth


def random_power(e_min, e_max, N):
    '''
    Returns random numbers from power law with given index and energy range (e_min, e_max in TeV)
    '''
    index = 2.48
    u = np.random.uniform(0, 1, N)
    # inverse transform sampling
    return ((e_max**(1 - index) - e_min**(1 - index)) * u + e_min**(1-index))**(1 / (1 - index))


@u.quantity_input(T=u.s, e_min=u.TeV, e_max=u.TeV, simulation_area=u.m**2)
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

    # Spectrum used in ctools
    index = 2.48
    C = u.Quantity(5.7e-6,  1 / (u.m**2 * u.s * u.TeV))
    E_0 = u.Quantity(0.3, u.TeV)

    return int((simulation_area.to(u.cm**2) * C * T * E_0**(index) / (1 - index)) * ((e_max)**(1 - index) - (e_min)**(1 - index)))


@u.quantity_input(T=u.s, e_min=u.TeV, e_max=u.TeV, simulation_area=u.m**2)
def number_particles_trans(T, e_min, e_max, simulation_area, z):
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
    C = u.Quantity(5.7e-6,  1 / (u.m**2 * u.s * u.TeV))
    E_0 = u.Quantity(0.3, u.TeV)

    tau = OptDepth.readmodel(model='dominguez')

    def spectrum_crab_ebl(z):
        return lambda E: C * (E / E_0) ** (-index) * np.exp(-1.0 * tau.opt_depth(z, E.to_value(u.TeV)))

    E = np.linspace(e_min, e_max, 20000)

    I = trapz(spectrum_crab_ebl(z)(E), E)

    return int(simulation_area * T * I)
