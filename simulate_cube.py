import pandas as pd
import numpy as np
import click
import astropy.units as u

from collections import OrderedDict
from ctawave.toy_models_crab import simulate_steady_source_with_transient, simulate_steady_source
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
from scipy import signal

from IPython import embed


'''
Simulation (cubes) for a transient in the field of view of a steady source.
3 simulation parts per cube: steady source, steady source with transient, steady source. Every part contains 'num_slices' images/slices.

Returns:
astropy tables:
1. cube_table: 'n_transient' cubes with  size [3*num_slices, bins_, bins_]
    transient template
2. trans_table: 'n_transient' timeseries with size [3*num_slices],
    crab units of flare (brightness)
    start/end-slice of transient ('simulation truth'),
    transient template
'''


@click.command()
@click.option(
    'output_path',
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    'irf_path',
    '--irf_path',
    type=click.Path(dir_okay=True),
    help='Directory for CTA Instrument Response Function (prod3b)',
    default='/home/lena/Dokumente/CTA'
)
@click.option(
    '--n_transient',
    '-n',
    type=click.INT,
    help='Number of transients to simulate',
    default='10'
)
@click.option(
    '--transient_template_index',
    '-temp',
    type=click.INT,
    help='Transient template index: 0=pks, 1=hess, 2=gauss',
    default='2'
)
@click.option(
    '--time_per_slice',
    '-t',
    type=click.INT,
    help='Observation time per slice in seconds',
    default='30'
)
@click.option(
    '--num_slices',
    '-s',
    type=click.INT,
    help='Number of slices per simulation. Size cube = 3*num_slices!!!!',
    default='20'
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
    n_transient,
    transient_template_index,
    time_per_slice,
    num_slices,
    bins_,
    cu_min=0.1,
    cu_max=7,
    duration_min=10,
    duration_max=100,
):
    '''
    Load CTA IRFs and transient template
    '''
    cta_perf_fits = fits.open('{}/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits'.format(irf_path))
    data_A_eff = cta_perf_fits['EFFECTIVE AREA']
    data_ang_res = cta_perf_fits['POINT SPREAD FUNCTION']
    data_bg_rate = cta_perf_fits['BACKGROUND']

    pks_data = np.loadtxt('data/PKS2155-flare06.dat', unpack=True)
    hess_data = np.loadtxt('data/LAT-GRB130427.dat', unpack=True)

    gauss = signal.gaussian(num_slices, std=1)
    transient_templates = [pks_data[1], hess_data[1], gauss]

    a_eff_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_A_eff.data['ENERG_LO'][0] + data_A_eff.data['ENERG_HI'][0])/2, "A_eff": data_A_eff.data['EFFAREA'][0][0]}))
    ang_res_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_ang_res.data['ENERG_LO'][0] + data_ang_res.data['ENERG_HI'][0])/2, "Ang_Res": data_ang_res.data['SIGMA_1'][0][0]}))


    ### used as range for random transient positions
    fov_min = 0 * u.deg
    fov_max = 12 * u.deg

    '''
    Start simulation
    '''
    slices = []
    trans_scales = []
    list_cu_flare = []
    for i in tqdm(range(n_transient)):

        '''
        Simulate slices containing one steady source and no transient
        '''
        slices_steady_source = simulate_steady_source(
                    x_pos=6*u.deg,
                    y_pos=6*u.deg,
                    df_A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    df_Ang_Res=ang_res_cta_south,
                    num_slices=num_slices,
                    time_per_slice=time_per_slice * u.s,
                    bins=[bins_, bins_],
                    )

        slices = np.append(slices, slices_steady_source)
        trans_scales = np.append(trans_scales, np.zeros(num_slices))

        '''
        Simulate slices with steady source and transient
        '''
        cu_flare = (cu_max - cu_min) * np.random.random() + cu_min
        list_cu_flare.append(cu_flare)

        slices_transient, trans_scale = simulate_steady_source_with_transient(
                    x_pos_steady_source=6*u.deg,
                    y_pos_steady_source=6*u.deg,
                    # x_pos_transient=np.random.randint(fov_min/u.deg, fov_max/u.deg)*u.deg,
                    # y_pos_transient=np.random.randint(fov_min/u.deg, fov_max/u.deg)*u.deg,
                    x_pos_transient=2*u.deg,
                    y_pos_transient=2*u.deg,
                    df_A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    df_Ang_Res=ang_res_cta_south,
                    cu_flare=cu_flare,
                    transient_template=transient_templates[transient_template_index],
                    # num_slices=np.random.randint(duration_min, duration_max),
                    num_slices=num_slices,
                    time_per_slice=time_per_slice * u.s,
                    bins=[bins_, bins_],
                    )
        trans_scales = np.append(trans_scales, trans_scale)
        slices = np.append(slices, slices_transient)

        '''
        Simulate slices containing one steady source and no transient
        '''
        slices_steady_source = simulate_steady_source(
                    x_pos=6*u.deg,
                    y_pos=6*u.deg,
                    df_A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    df_Ang_Res=ang_res_cta_south,
                    num_slices=num_slices,
                    time_per_slice=time_per_slice * u.s,
                    bins=[bins_, bins_],
                    )
        slices = np.append(slices, slices_steady_source)
        trans_scales = np.append(trans_scales, np.zeros(num_slices))

    '''
    Write and save astropy tables for simulated cubes and transients.
    '''

    list_cubes = slices.reshape([-1, 3*num_slices, bins_, bins_])
    list_transients = trans_scales.reshape([-1, 3*num_slices])

    cube_table = Table()
    trans_table = Table()

    trans_table['timeseries'] = list_transients
    trans_table['cu_flare'] = list_cu_flare
    trans_table['template'] = transient_template_index

    ### start, end slice for gaussian shape, arbitrary!!
    trans_table['start_flare'] = 7 + num_slices
    trans_table['end_flare'] = 12 + num_slices

    cube_table['cube'] = list_cubes
    cube_table['template'] = transient_template_index
    cube_table['num_slices'] = 3*num_slices
    cube_table['num_flare_slices'] = num_slices

    cube_table.meta['n_transient'] = n_transient
    cube_table.meta['num_slices'] = 3*num_slices
    cube_table.meta['template'] = transient_template_index
    cube_table.meta['time_per_slice'] = time_per_slice
    cube_table.meta['bins'] = bins_

    trans_table.meta['n_transient'] = n_transient
    trans_table.meta['num_slices'] = 3*num_slices
    trans_table.meta['template'] = transient_template_index
    trans_table.meta['time_per_slice'] = time_per_slice
    trans_table.meta['bins'] = bins_

    cube_table.write('{}/n{}_s{}_t{}_cube.hdf5'.format(output_path, n_transient, 3*num_slices, transient_template_index), path='data', overwrite=True)
    trans_table.write('{}/n{}_s{}_t{}_trans.hdf5'.format(output_path, n_transient, 3*num_slices, transient_template_index), path='data', overwrite=True)


if __name__ == '__main__':
    main()
