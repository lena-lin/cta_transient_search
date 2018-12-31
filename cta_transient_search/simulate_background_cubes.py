import click
import astropy.units as u
import numpy as np
from collections import OrderedDict
from simulation.simulate_skymap import simulate_steady_source
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm


@click.command()
@click.option(
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--irf_path',
    type=click.Path(dir_okay=True),
    help='Directory for CTA Instrument Response Function (prod3b)',
    default='/home/lena/Dokumente/CTA'
)
@click.option(
    '--n_cubes',
    '-n',
    type=click.INT,
    help='Number of cubes to simulate',
    default='10'
)
@click.option(
    '--num_slices',
    '-s',
    type=click.INT,
    help='Number of slices per cube',
    default='60'
)
@click.option(
    '--time_per_slice',
    '-t',
    type=click.INT,
    help='Time per slice',
    default='30'
)
@click.option(
    '--bins_',
    '-b',
    type=click.INT,
    help='Number of ra/dec bins for simulated cube',
    default='80'
)
def main(
    output_path,
    irf_path,
    n_cubes,
    num_slices,
    time_per_slice,
    bins_,
):
    cta_perf_fits = fits.open('{}/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits'.format(irf_path))
    data_A_eff = cta_perf_fits['EFFECTIVE AREA']
    data_ang_res = cta_perf_fits['POINT SPREAD FUNCTION']
    data_bg_rate = cta_perf_fits['BACKGROUND']

    a_eff_cta_south = OrderedDict(
                            {
                                "E_TeV": (data_A_eff.data['ENERG_LO'][0] + data_A_eff.data['ENERG_HI'][0])/2,
                                "A_eff": data_A_eff.data['EFFAREA'][0]
                            }
                        )

    psf_cta_south = OrderedDict(
                                {
                                    "E_TeV": (data_ang_res.data['ENERG_LO'][0] + data_ang_res.data['ENERG_HI'][0])/2,
                                    "psf_sigma": data_ang_res.data['SIGMA_1'][0]
                                }
                            )

    slices = []
    for i in tqdm(range(n_cubes)):

        slices_steady_source = simulate_steady_source(
                    A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    psf=psf_cta_south,
                    num_slices=num_slices,
                    time_per_slice=time_per_slice * u.s,
                    bins=[bins_, bins_],
                    )

        slices.append(slices_steady_source)

    bg_table = Table()
    bg_table.meta['n_cubes'] = n_cubes
    bg_table.meta['num_slices'] = num_slices
    bg_table.meta['time_per_slice'] = time_per_slice
    bg_table.meta['bins'] = bins_
    bg_table.meta['fov'] = 12 * u.deg
    bg_table.meta['steady_source'] = 'Crab'

    bg_table['cube'] = slices
    bg_table.write('{}/n{}_s{}_bg.hdf5'.format(output_path, n_cubes, num_slices), path='data', overwrite=True)


if __name__ == '__main__':
    main()
