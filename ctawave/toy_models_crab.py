import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import spectrum
import performance


def simulate_steady_source_with_transient(
            A_eff,
            fits_bg_rate,
            psf,
            cu_flare,
            pos_ra,
            pos_dec,
            transient_template,
            num_slices=100,
            time_per_slice=30 * u.s,
            bins=[80, 80],
            E_min=0.1 * u.TeV,
            E_max=100 * u.TeV,
            fov=12 * u.deg,
            ):

    cta_radius = 800 * u.m
    sim_area = 2*np.pi*(2*cta_radius)**2
    crab_coord = SkyCoord.from_name('Crab')

    N_steady_source = spectrum.number_particles_crab(time_per_slice, E_min, E_max, sim_area)
    N_background_cta = performance.integrate_background(fits_bg_rate, time_per_slice)

    flare_interp = np.interp(range(num_slices), np.linspace(0, num_slices, len(transient_template)), transient_template)
    transient_scale = (flare_interp/flare_interp.max() * N_steady_source*cu_flare).astype(int)

    # valid_transient_position = False
    # while not valid_transient_position:
    #     ra_transient = np.random.randint(
    #                                     crab_coord.ra.deg - fov.value / 2 + fov.value / 10,
    #                                     crab_coord.ra.deg + fov.value / 2 - fov.value / 10
    #                                     )
    #     dec_transient = np.random.randint(
    #                                     crab_coord.dec.deg - fov.value / 2 + fov.value / 10,
    #                                     crab_coord.dec.deg + fov.value / 2 - fov.value / 10
    #                                     )
    #     theta = np.sqrt((crab_coord.ra.deg - ra_transient)**2 + (crab_coord.dec.deg - dec_transient)**2)
    #     if theta > 1 and theta < 6:
    #         valid_transient_position = True

    ra_transient = crab_coord.ra.deg - fov.value / pos_ra
    dec_transient = crab_coord.dec.deg - fov.value / pos_dec
    theta = np.sqrt((crab_coord.ra.deg - ra_transient)**2 + (crab_coord.dec.deg - dec_transient)**2)

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
                                                fits_bg_rate,
                                                time_per_slice,
                                                N_background_cta,
                                                fov,
                                                crab_coord,
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

        else:
            RA_tr, DEC_tr = [], []

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

    return np.array(slices), transient_scale, ra_transient, dec_transient


def simulate_steady_source(
            A_eff,
            fits_bg_rate,
            psf,
            num_slices,
            time_per_slice=30 * u.s,
            bins=[80, 80],
            E_min=0.1 * u.TeV,
            E_max=100 * u.TeV,
            fov=12 * u.deg,
        ):
    cta_radius = 800 * u.m
    sim_area = 2*np.pi*(2*cta_radius)**2
    crab_coord = SkyCoord.from_name('Crab')

    N_steady_source = spectrum.number_particles_crab(time_per_slice, E_min, E_max, sim_area)
    N_background_cta = performance.integrate_background(fits_bg_rate, time_per_slice)
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
                                                fits_bg_rate,
                                                time_per_slice,
                                                N_background_cta,
                                                fov,
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
