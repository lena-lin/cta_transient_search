import pandas as pd
import numpy as np
import click
import astropy.units as u

from collections import OrderedDict
from ctawave.toy_models_crab import simulate_steady_source
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

from IPython import embed


@click.command()
@click.argument('out_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '--n_cubes',
    '-n',
    type=click.INT,
    help='Number of cubes to simulate',
    default='200'
)
@click.option(
    '--n_slices',
    '-s',
    type=click.INT,
    help='Number of slices per cube',
    default='10'
)
@click.option(
    '--time_per_slice',
    '-t',
    type=click.INT,
    help='Time per slice',
    default='30'
)


def main(
    out_file,
    n_cubes,
    n_slices,
    time_per_slice
    ):
    cta_perf_fits = fits.open('/home/lena/Software/ctools-1.3.0/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits')
    data_A_eff = cta_perf_fits['EFFECTIVE AREA']
    data_ang_res = cta_perf_fits['POINT SPREAD FUNCTION']
    data_bg_rate = cta_perf_fits['BACKGROUND']

    a_eff_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_A_eff.data['ENERG_LO'][0] + data_A_eff.data['ENERG_HI'][0])/2, "A_eff": data_A_eff.data['EFFAREA'][0][0]}))
    ang_res_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_ang_res.data['ENERG_LO'][0] + data_ang_res.data['ENERG_HI'][0])/2, "Ang_Res": data_ang_res.data['SIGMA_1'][0][0]}))
    print('start')
    slices = []
    for i in tqdm(range(n_cubes)):
        slices_steady_source = simulate_steady_source(
                    x_pos=6*u.deg,
                    y_pos=6*u.deg,
                    df_A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    df_Ang_Res=ang_res_cta_south,
                    num_slices=n_slices
                    )
        slices.append(slices_steady_source)

    embed()
    bg_table = Table()
    bg_table['bg_cubes'] = slices
    bg_table.write('build/n200_s60_bg.hdf5', path='data')




if __name__ == '__main__':
    print('main')
    main()
