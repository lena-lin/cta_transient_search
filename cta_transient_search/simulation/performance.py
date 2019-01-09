import numpy as np
from .spectrum import random_power
import astropy.units as u
from scipy import integrate
import pandas as pd
from scipy.interpolate import interp1d, splrep, splev
# from IPython import embed


def interp_ang_res(event_E_TeV, bin_bound_E_TeV, psf_sigma):
    '''
    Interpolates the angular resolution for given energy E_Tev
    '''
    return np.interp(event_E_TeV, bin_bound_E_TeV, psf_sigma)


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
    Returns array of events (list of energies) from a source with power law distribution
    folded with effective area of the telescope

    Parameters:
    index: index of the power law spectrum
    e_min: lower energy bound
    e_max: upper energy bound
    N: number of events coming from source
    A_eff: DataFrame containing values for effective area
    '''

    a_eff_func = interp1d(A_eff['E_TeV'], A_eff['A_eff'][theta])

    simulation_area = simulation_area.to(u.m**2)
    events = random_power(e_min, e_max, N)
    folded_events = []
    if len(events) > 0:
        for e in events:
            a_eff_event = a_eff_func(e)
            ulimite = (a_eff_event * u.m**2) / (simulation_area)

            if(ulimite.value >= np.random.uniform(0, 1)):
                folded_events.append(e/u.TeV)

    return folded_events


def sample_positions_steady_source(ra_pos, dec_pos, ang_res):
    '''
    Sample position for every particle with given mean position and angular resolution as
    mean and covariance for normal distribution
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


# def sample_positions_background_random(irf_bg, time_per_slice, N, fov, crab_coord):
#     '''
#     Sample positions for given number of background events from background distribution
#     '''
#     ra_min = crab_coord.ra.deg - fov.value/2
#     ra_max = crab_coord.ra.deg + fov.value/2
#     dec_min = crab_coord.dec.deg - fov.value/2
#     dec_max = crab_coord.dec.deg + fov.value/2
#
#     bg = irf_bg.data['BGD'][0]
#     bg_dist_for_slice = bg.sum(axis=0) * time_per_slice.value
#
#     bins_ra = bg_dist_for_slice.shape[0]
#     bins_dec = bg_dist_for_slice.shape[1]
#
#     cumsum_cols = np.cumsum(bg_dist_for_slice, axis=0)
#     cumsum_last_row = np.cumsum(cumsum_cols[-1])
#     y = np.linspace(dec_min, dec_max, bins_dec + 1)
#
#     inverse_response_curves_y = []
#     for col in cumsum_cols.T:
#         col = np.hstack([0, col])
#         inverse_response_curves_y.append(interpolate.interp1d(col, y))
#
#     y = np.linspace(ra_min, ra_max, bins_ra + 1)
#     inverse_response_x = interpolate.interp1d(np.hstack([0, cumsum_last_row]), y)
#
#     RA_bg = []
#     DEC_bg = []
#     for i in range(int(N.value)):
#         urx = np.random.uniform(0, cumsum_last_row.max())
#         ra = inverse_response_x(urx)
#         ix = int((ra - ra_min)/(ra_max - ra_min) * bins_ra)
#         ury = np.random.uniform(0, cumsum_cols.T[ix].max())
#         dec = inverse_response_curves_y[ix](ury)
#         RA_bg.append(ra)
#         DEC_bg.append(dec)
#
#     return RA_bg, DEC_bg

def sample_positions_background_random(N, inv_F, cumsum_min, crab_coord):
    u = np.random.uniform(cumsum_min, 1, N)
    r_sampled = inv_F(u)

    phi = np.random.uniform(0, 2 * np.pi, len(u))

    ra_bg = r_sampled * np.cos(phi) + crab_coord.ra.deg
    dec_bg = r_sampled * np.sin(phi) + crab_coord.dec.deg

    return ra_bg, dec_bg


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
    return int((integral.sum() * obs_time).value)


def inverse_response_background(irf):
    bg = irf['BACKGROUND'].data['BGD'][0].sum(axis=0)
    bg = bg.flatten()

    x = (irf['BACKGROUND'].data['DETX_HI'] + irf['BACKGROUND'].data['DETX_LO'])/2
    y = (irf['BACKGROUND'].data['DETY_HI'] + irf['BACKGROUND'].data['DETY_LO'])/2

    X, Y = np.meshgrid(x, y)
    r = np.sqrt(X ** 2 + Y ** 2)
    r = r.flatten()
    sort_idx = np.argsort(r)
    r = r[sort_idx]
    bg = bg[sort_idx]

    bg_cut = bg[:np.where(bg == 0)[0].min()]
    r_cut = r[:np.where(bg == 0)[0].min()]

    df_bg = pd.DataFrame(data={'r': r_cut, 'bg': bg_cut})
    df_bg_grouped = df_bg.groupby(['r']).mean()

    r_unique = df_bg_grouped.index.values
    bg_mean = df_bg_grouped['bg'].values

    bg_mean = np.append(bg_mean[0], bg_mean)
    r_unique = np.append(0, r_unique)

    r_even = np.linspace(r_unique.min(), r_unique.max(), 100)
    spline_f = splrep(r_unique, r_unique * bg_mean, xb=0, xe=6, k=3)
    cumsum = splev(r_even, spline_f).cumsum()/splev(r_even, spline_f).cumsum().max()
    inv_F = interp1d(cumsum, r_even)

    return inv_F, cumsum.min()
