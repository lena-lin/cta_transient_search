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

from LC_forms import Simple_Gaussian, Small_Gaussian, Exponential
from random import randint


'''
Simulation (cubes) for a transient in the field of view of a steady source.
3 simulation parts per cube: steady source, steady source with transient, steady source. Every part contains 'num_slices' images/slices.
Input:
See Click.options()

Returns:
astropy tables:
1. cube_table: 'n_transient' cubes with  size [3*num_slices, bins_, bins_]
    transient template
2. trans_table: 'n_transient' timeseries with size [3*num_slices],
    crab units of flare (brightness)
    start/end-slice of transient ('simulation truth'),
    transient template
'''

def pull_template_index(x):
    x = randint(2,5)
    return x



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
    default=''
)
@click.option(
    '--n_transient',
    '-n',
    type=click.INT,
    help='Number of transients to simulate',
    default='10'
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
@click.option(
    '--percentage_transient_noise',
    '-p',
    type=click.INT,
    help='strength of randomized noise on transient template in percent',
    default='4'
)
def main(
    output_path,
    irf_path,
    n_transient,
    time_per_slice,
    num_slices,
    bins_,
    percentage_transient_noise, #**
    cu_min=1, #0.1
    cu_max=7,#7
    duration_min=10,
    duration_max=100,
):
    '''
    Load CTA IRFs and transient template
    '''
    cta_perf_fits = fits.open('/home/jana/Schreibtisch/Projekt Master/Masterarbeit/Templates/Sensitivity/IRFs/South_z20_average_100s/irf_file.fits')
    data_A_eff = cta_perf_fits['EFFECTIVE AREA']
    data_ang_res = cta_perf_fits['POINT SPREAD FUNCTION']
    data_bg_rate = cta_perf_fits['BACKGROUND']
# First Choise of used templates to interpolate
    pks_data = np.loadtxt('data/PKS2155-flare06.dat', unpack=True)
    hess_data = np.loadtxt('data/LAT-GRB130427.dat', unpack=True)
# simple gaussian, std= 1
    gauss = signal.gaussian(num_slices, std=1)
# new Templates after fitting gaussian + exponential to data
    simple = Simple_Gaussian(num_slices,percentage_transient_noise) # 4% noise
    small = Small_Gaussian(num_slices,percentage_transient_noise)
    exponential = Exponential(num_slices,percentage_transient_noise)

    transient_templates = [pks_data[1], hess_data[1], gauss,simple,small,exponential]  # indices 0 to 5
# Choose start of transient dependent on template
    transient_start_slices = np.array([20,20,num_slices/2.0-3, num_slices/2.0-5, num_slices/2.0-1, num_slices/2.0-1.0/3.0*num_slices])

    a_eff_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_A_eff.data['ENERG_LO'][0] + data_A_eff.data['ENERG_HI'][0])/2, "A_eff": data_A_eff.data['EFFAREA'][0][0]}))
    ang_res_cta_south = pd.DataFrame(OrderedDict({"E_TeV": (data_ang_res.data['ENERG_LO'][0] + data_ang_res.data['ENERG_HI'][0])/2, "Ang_Res": data_ang_res.data['SIGMA_1'][0][0]}))


    ### used as range for random transient positions
    fov_min = 0 * u.deg
    fov_max = 12 * u.deg

    '''
    Pull transient template indices
    '''
    transient_template_index = list(map(pull_template_index,np.zeros(n_transient)))  #***

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
        list_cu_flare.append(cu_flare) ## nicht vrest  mitnehmen?

        print('Start simulating transient with template_index', transient_template_index[i]) # **
        slices_transient, trans_scale = simulate_steady_source_with_transient(
                    x_pos_steady_source=6*u.deg,
                    y_pos_steady_source=6*u.deg,
                    # x_pos_transient=np.random.randint(fov_min/u.deg, fov_max/u.deg)*u.deg,
                    # y_pos_transient=np.random.randint(fov_min/u.deg, fov_max/u.deg)*u.deg,
                    x_pos_transient=5*u.deg,
                    y_pos_transient=5*u.deg,
                    df_A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    df_Ang_Res=ang_res_cta_south,
                    cu_flare=cu_flare,
                    transient_template=transient_templates[transient_template_index[i]],    # ***
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
    #print('Length Slices: ', len(slices))
    #print('Length Cu_Flare List', len(list_cu_flare))
    #print('Length Transient scales', len(trans_scales))

    list_cubes = slices.reshape([-1, 3*num_slices, bins_, bins_])
    list_transients = trans_scales.reshape([-1, 3*num_slices])

    cube_table = Table()
    trans_table = Table()

    trans_table['timeseries'] = list_transients
    trans_table['cu_flare'] = list_cu_flare
    trans_table['template'] = transient_template_index

    ### start slice for templates, dependent on template index + num_slices
    ### end slice arbitrary!! Not used
    trans_table['start_flare'] = transient_start_slices[transient_template_index] + num_slices  ## default 7
    trans_table['end_flare'] = 12 + num_slices

    cube_table['cube'] = list_cubes
    cube_table['template'] = transient_template_index  ## Bisher eine Zahl ?
    cube_table['num_slices'] = 3*num_slices
    cube_table['num_flare_slices'] = num_slices

    cube_table.meta['n_transient'] = n_transient
    cube_table.meta['percentage_transient_noise'] = percentage_transient_noise
    cube_table.meta['num_slices'] = 3*num_slices
    cube_table.meta['template'] = transient_template_index
    cube_table.meta['time_per_slice'] = time_per_slice
    cube_table.meta['bins'] = bins_

    trans_table.meta['n_transient'] = n_transient
    trans_table.meta['percentage_transient_noise'] = percentage_transient_noise
    trans_table.meta['num_slices'] = 3*num_slices
    trans_table.meta['template'] = transient_template_index
    trans_table.meta['time_per_slice'] = time_per_slice
    trans_table.meta['bins'] = bins_

    cube_table.write('{}/n{}_s{}_p{}_cube.hdf5'.format(output_path, n_transient, 3*num_slices,percentage_transient_noise), path='data', overwrite=True)
    trans_table.write('{}/n{}_s{}_p{}_trans.hdf5'.format(output_path, n_transient, 3*num_slices,percentage_transient_noise), path='data', overwrite=True)


if __name__ == '__main__':
    main()
