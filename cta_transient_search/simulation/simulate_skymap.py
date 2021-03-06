import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from .spectrum import number_particles_crab, number_particles_trans
from simulation import performance


def simulate_steady_source_with_transient(
            A_eff,
            N_background_cta,
            inv_F,
            cumsum_min,
            psf,
            cu_flare,
            pos_ra,
            pos_dec,
            pos_random,
            transient_template,
            z_trans,
            num_slices,
            time_per_slice=10 * u.s,
            bins=[80, 80],
            E_min=0.1 * u.TeV,
            E_max=100 * u.TeV,
            fov=8 * u.deg,
            ):

    cta_radius = 800 * u.m
    sim_area = 2*np.pi*(2*cta_radius)**2
    #crab_coord = SkyCoord.from_name('Crab')
    crab_coord = SkyCoord('05 34 31.97 +22 00 52.1', unit=(u.hourangle, u.deg))
    N_steady_source = number_particles_crab(time_per_slice, E_min, E_max, sim_area)
    N_transient_ebl = number_particles_trans(time_per_slice, E_min, E_max, sim_area, z_trans)

    flare_interp = np.interp(range(num_slices), np.linspace(0, num_slices, len(transient_template)), transient_template)
    transient_scale = (flare_interp/flare_interp.max() * N_transient_ebl*cu_flare).astype(int)

    if pos_random == True:
        valid_transient_position = False
        while not valid_transient_position:
            ra_transient = np.random.randint(
                                         crab_coord.ra.deg - fov.value / 2 + fov.value / 10,
                                         crab_coord.ra.deg + fov.value / 2 - fov.value / 10
                                         )
            dec_transient = np.random.randint(
                                         crab_coord.dec.deg - fov.value / 2 + fov.value / 10,
                                         crab_coord.dec.deg + fov.value / 2 - fov.value / 10
                                         )
            theta = np.sqrt((crab_coord.ra.deg - ra_transient)**2 + (crab_coord.dec.deg - dec_transient)**2)
            if theta > 1 and theta < 4:
                valid_transient_position = True

    if pos_random == False:
        ra_transient = crab_coord.ra.deg - fov.value / pos_ra
        dec_transient = crab_coord.dec.deg - fov.value / pos_dec
        theta = np.sqrt((crab_coord.ra.deg - ra_transient)**2 + (crab_coord.dec.deg - dec_transient)**2)

    slices = []
    num_transient_events = []
    for i in range(num_slices):
        folded_events_crab = performance.response(
                                                time_per_slice,
                                                N_steady_source,
                                                E_min, E_max,
                                                A_eff,
                                                sim_area,
                                                theta=0,
                                            )
        ang_res_steady_source = performance.interp_ang_res(
                                                folded_events_crab,
                                                psf['E_TeV'],
                                                psf['psf_sigma'][0]
                                            )
        RA_crab, DEC_crab = performance.sample_positions_steady_source(
                                                crab_coord.ra.deg,
                                                crab_coord.dec.deg,
                                                ang_res_steady_source,
                                            )
        RA_bg, DEC_bg = performance.sample_positions_background_random(
                                                N_background_cta,
                                                inv_F,
                                                cumsum_min,
                                                crab_coord
                                            )

        if transient_scale[i] > 0:
            folded_events_transient = performance.response(
                                                    time_per_slice,
                                                    transient_scale[i],
                                                    E_min,
                                                    E_max,
                                                    A_eff,
                                                    sim_area,
                                                    theta=int(theta)
                                                )
            ang_res_transinet = performance.interp_ang_res(
                                                    folded_events_transient,
                                                    psf['E_TeV'],
                                                    psf['psf_sigma'][int(theta)]
                                                )
            RA_tr, DEC_tr = performance.sample_positions_steady_source(
                                                    ra_transient,
                                                    dec_transient,
                                                    ang_res_transinet,
                                                )
            num_transient_events.append(len(folded_events_transient))

        else:
            RA_tr, DEC_tr = [], []
            num_transient_events.append(0)

        RA = np.concatenate([RA_bg, RA_tr, RA_crab])
        DEC = np.concatenate([DEC_bg, DEC_tr, DEC_crab])

        slices.append(
                np.histogram2d(
                    RA,
                    DEC,
                    range=[
                        [crab_coord.ra.deg - fov.value / 2, crab_coord.ra.deg + fov.value / 2],
                        [crab_coord.dec.deg - fov.value / 2, crab_coord.dec.deg + fov.value / 2]
                    ],
                    bins=bins,
                )[0]
        )



    return np.array(slices), num_transient_events, ra_transient, dec_transient


def simulate_steady_source(
            A_eff,
            N_background_cta,
            inv_F,
            cumsum_min,
            psf,
            num_slices,
            time_per_slice=10 * u.s,
            bins=[80, 80],
            E_min=0.1 * u.TeV,
            E_max=100 * u.TeV,
            fov=8 * u.deg,
        ):
    cta_radius = 800 * u.m
    sim_area = 2*np.pi*(2*cta_radius)**2
    #crab_coord = SkyCoord.from_name('Crab')  ## Astropy exception handling in 'from_name', BUGG :(
    crab_coord = SkyCoord('05 34 31.97 +22 00 52.1', unit=(u.hourangle, u.deg))

    N_steady_source = number_particles_crab(time_per_slice, E_min, E_max, sim_area)
    n_events = 0

    slices = []
    for i in range(num_slices):
        folded_events_crab = performance.response(
                                            time_per_slice,
                                            N_steady_source,
                                            E_min, E_max,
                                            A_eff,
                                            sim_area,
                                            theta=0,
                                        )
        ang_res_steady_source = performance.interp_ang_res(
                                                folded_events_crab,
                                                psf['E_TeV'],
                                                psf['psf_sigma'][0]
                                        )

        RA_crab, DEC_crab = performance.sample_positions_steady_source(
                                            crab_coord.ra.deg,
                                            crab_coord.dec.deg,
                                            ang_res_steady_source,
                                        )
        RA_bg, DEC_bg = performance.sample_positions_background_random(
                                                N_background_cta,
                                                inv_F,
                                                cumsum_min,
                                                crab_coord
                                            )
        RA = np.concatenate([RA_bg, RA_crab])
        DEC = np.concatenate([DEC_bg, DEC_crab])

        slices.append(
                np.histogram2d(
                    RA,
                    DEC,
                    range=[
                        [crab_coord.ra.deg - fov.value / 2, crab_coord.ra.deg + fov.value / 2],
                        [crab_coord.dec.deg - fov.value / 2, crab_coord.dec.deg + fov.value / 2]
                    ],
                    bins=bins,
                )[0]
        )
        n_events += 1/float(num_slices)*len(folded_events_crab)

    return np.array(slices)
