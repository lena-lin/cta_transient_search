import numpy as np
import spectrum
import astropy.units as u
from scipy import integrate
from IPython import embed


def interp_ang_res(event_E_TeV, bin_bound_E_TeV, psf_sigma):
    '''
    Interpolates the angular resolution for given energy E_Tev
    '''
    return np.interp(event_E_TeV, bin_bound_E_TeV, psf_sigma)

#TODO check calc_a_eff_factor
def calc_a_eff_factor(df_A_eff, simulation_area):
    '''
    Returns integrated effective area divided by full impact area
    '''

    simulation_area = simulation_area.to(u.m**2)
    integrate_a_eff = integrate.simps(y=df_A_eff.A_eff, x=df_A_eff.E_TeV) * u.m**2
    integrate_a_impact = simulation_area*(df_A_eff.E_TeV.max() - df_A_eff.E_TeV.min())
    return integrate_a_eff/integrate_a_impact


def interp_eff_area(event_E_TeV, bin_bound_E_TeV, a_eff):
    '''
    Interpolates effective area for given energy
    '''
    return np.interp(event_E_TeV / u.TeV, bin_bound_E_TeV, a_eff)


def response(T, N, e_min, e_max, A_eff, simulation_area, theta):
    '''
    Returns array of events (list of energies) from a source with power law distribution folded with effective area of the telescope

    Parameters:
    index: index of the power law spectrum
    e_min: lower energy bound
    e_max: upper energy bound
    N: number of events coming from source
    A_eff: DataFrame containing values for effective area
    '''

    simulation_area = simulation_area.to(u.m**2)
    events = spectrum.random_power(e_min, e_max, N)
    folded_events = []
    if len(events) > 0:
        for e in events:
            a_eff_event = interp_eff_area(e, A_eff['E_TeV'], A_eff['A_eff'][theta])
            ulimite = (a_eff_event * u.m**2) / (simulation_area)

            if(ulimite.value >= np.random.uniform(0, 1)):
                folded_events.append(e/u.TeV)

    return folded_events


def sample_positions_steady_source(ra_pos, dec_pos, ang_res):
    '''
    Sample position for every particle with given mean position and angular resolution as mean and covariance for normal distribution
    '''
    mean = [ra_pos, dec_pos]
    RA = []
    DEC = []
    for r in ang_res:
        cov = [[r**2, 0], [0, r**2]]
        x, y = np.random.multivariate_normal(mean, cov).T
        RA.append(x)
        DEC.append(y)
    return RA, DEC


def sample_positions_background_random(fov, source_coordinates, N):
    '''
    Sample positions for given number of background events from normal distribution
    '''
    RA_bg = np.random.uniform(
                        source_coordinates.ra.value - fov.value / 2,
                        source_coordinates.ra.value + fov.value / 2,
                        N,
                    )
    DEC_bg = np.random.uniform(
                        source_coordinates.dec.value - fov.value / 2,
                        source_coordinates.dec.value + fov.value / 2,
                        N,
                    )

    return RA_bg, DEC_bg


def integrate_background(bkg, obs_time):
    '''
    Integrate background events from IRF.
    fits-file provides background rates in 1/MeV/s/sr, therefore energy units must be adapted!
    There are different IRFs for different observation times and zenith angles. Here the IRF for
    20 deg and 100 s is used.

    Parameters:
    bkg: irf-file containing background cube
    obs_time: observation time per slice
    '''

    delta_energy = (bkg.data['ENERG_HI'][0] - bkg.data['ENERG_LO'][0]) * u.TeV
    delta_x = (bkg.data['DETX_HI'][0] - bkg.data['DETX_LO'][0]) * u.deg
    delta_y = (bkg.data['DETY_HI'][0] - bkg.data['DETY_LO'][0]) * u.deg
    delta_energy, delta_y, delta_x = np.meshgrid(delta_energy, delta_y, delta_x, indexing='ij')
    bin_volume = delta_energy.to(u.MeV) * (delta_y * delta_x).to(u.sr)
    bg = bkg.data['BGD'][0] * (1/(u.MeV * u.s * u.sr))
    integral = bg * bin_volume
    return integral.sum() * obs_time
