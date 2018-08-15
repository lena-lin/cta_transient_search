import numpy as np
import click
import astropy.units as u

from collections import OrderedDict
from ctawave.toy_models_crab import simulate_steady_source_with_transient, simulate_steady_source
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

# from IPython import embed

#from LC_forms import broad_gaussian, narrow_gaussian, deltapeak_exponential  # Old version
from New_LC_forms import simulate_Gaussians,simulate_Exponential # New test


'''
Simulation (cubes) for a transient in the field of view of a steady source.
3 simulation parts per cube: steady source, steady source with transient, steady source.
Every part contains 'num_slices' images/slices.

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
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--irf_path',
    '-f',
    type=click.Path(dir_okay=True),
    help='Directory for CTA Instrument Response Function (prod3b)',
    #default='/home/lena/Dokumente/CTA'
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
    '-i',
    type=click.INT,
    help='Transient template index:  0=broad gauss, 1=narrow gauss, 2=deltapeak + exponential fall',
    default='2'
)
@click.option(
    '--random_transient_template',
    '-r',
    help='Samples random template for transient shape, if True',
    is_flag=True
)
@click.option(
    '--noise',
    type=click.INT,
    help='Noise for transient template in %',
    default=0
)
@click.option(
    '--time_per_slice',
    '-t',
    type=click.INT,
    help='Observation time per slice in seconds',
    default='10'
)
@click.option(
    '--cu_min',
    type=click.FLOAT,
    help='Minimum transient brightness in crab units',
    default=8
)
@click.option(
    '--cu_max',
    type=click.FLOAT,
    help='Maximum transient brightness in crab units',
    default=10
)
@click.option(
    '--trans_pos_ra',
    '-x',
    type=click.INT,
    help='Ra Transient (Ration FoV!!): crab_coord.ra.deg - fov.value / X',
    default=8
)
@click.option(
    '--trans_pos_dec',
    '-y',
    type=click.INT,
    help='Dec Transient (Ration FoV!!): crab_coord.dec.deg - fov.value / Y',
    default=8
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
    noise,
    random_transient_template,
    time_per_slice,
    num_slices,
    bins_,
    cu_min,
    cu_max,
    trans_pos_ra,
    trans_pos_dec
):
    '''
    Load CTA IRFs and transient template
    '''
    cta_perf_fits = fits.open('{}/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits'.format(irf_path))
    data_A_eff = cta_perf_fits['EFFECTIVE AREA']
    data_ang_res = cta_perf_fits['POINT SPREAD FUNCTION']
    data_bg_rate = cta_perf_fits['BACKGROUND']
# First Choice of used templates to interpolate : Problem of defining the truth
    pks_data = np.loadtxt('data/PKS2155-flare06.dat', unpack=True)
    hess_data = np.loadtxt('data/LAT-GRB130427.dat', unpack=True)

# new Templates after fitting gaussian + exponential to data
    simple, true_start_simple = simulate_Gaussians(1.8348, 16.0364, num_slices, time_per_slice)
    small, true_start_small = simulate_Gaussians(0.45, 2.18, num_slices, time_per_slice)
    exponential, true_start_exponential = simulate_Exponential(3, 6, 0, 2, num_slices, time_per_slice)

    transient_templates = [simple, small, exponential]

# Choose start of transient dependent on template
    transient_start_slices = np.array([
                                        true_start_simple, true_start_small, true_start_exponential
                                        ])

    a_eff_cta_south = OrderedDict({
                                "E_TeV": (data_A_eff.data['ENERG_LO'][0] + data_A_eff.data['ENERG_HI'][0])/2,
                                "A_eff": data_A_eff.data['EFFAREA'][0]
                            })

    psf_cta_south = OrderedDict({
                                    "E_TeV": (data_ang_res.data['ENERG_LO'][0] + data_ang_res.data['ENERG_HI'][0])/2,
                                    "psf_sigma": data_ang_res.data['SIGMA_1'][0]
                                })

    '''
    Start simulation
    '''
    slices = []
    trans_scales = []
    list_cu_flare = []
    list_ra_transient = []
    list_dec_transient = []

    if random_transient_template is True:
        list_templates = np.random.randint(0, len(transient_templates), n_transient)
        transient_template_filename = 'random'
    else:
        try:
            list_templates = [transient_template_index] * n_transient
            transient_template_filename = transient_template_index
        except:
            print('Transient Template does not exist')

    for i in tqdm(range(n_transient)):

        '''
        Simulate slices containing one steady source and no transient
        '''
        slices_steady_source = simulate_steady_source(
                    A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    psf=psf_cta_south,
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

        slices_transient, trans_scale, ra_transient, dec_transient = simulate_steady_source_with_transient(
                    A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    psf=psf_cta_south,
                    cu_flare=cu_flare,
                    pos_ra=trans_pos_ra,
                    pos_dec=trans_pos_dec,
                    transient_template=transient_templates[list_templates[i]],
                    num_slices=num_slices,
                    time_per_slice=time_per_slice * u.s,
                    bins=[bins_, bins_],
                    )
        list_ra_transient.append(ra_transient)
        list_dec_transient.append(dec_transient)
        trans_scales = np.append(trans_scales, trans_scale)
        slices = np.append(slices, slices_transient)

        '''
        Simulate slices containing one steady source and no transient
        '''
        slices_steady_source = simulate_steady_source(
                    A_eff=a_eff_cta_south,
                    fits_bg_rate=data_bg_rate,
                    psf=psf_cta_south,
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
    list_transient_positions = np.dstack((list_ra_transient, list_dec_transient))[0]

    cube_table = Table()
    trans_table = Table()

    trans_table['timeseries'] = list_transients
    trans_table['cu_flare'] = list_cu_flare
    trans_table['template'] = list_templates
    trans_table['position'] = list_transient_positions
    # start slice for templates, dependent on template index + num_slices
    trans_table['start_flare'] = np.asanyarray([transient_start_slices[template] for template in list_templates]) + num_slices # add first empty cube, but not added in evaluation.py
    # end slice arbitrary!! Not used so far
    trans_table['end_flare'] = 12 + num_slices

    cube_table['cube'] = list_cubes
    cube_table['template'] = list_templates
    cube_table['num_flare_slices'] = num_slices

    cube_table.meta['n_transient'] = n_transient
    cube_table.meta['num_slices'] = 3*num_slices
    cube_table.meta['template'] = transient_template_filename
    cube_table.meta['time_per_slice'] = time_per_slice
    cube_table.meta['bins'] = bins_
    cube_table.meta['fov'] = 12 * u.deg
    cube_table.meta['steady_source'] = 'Crab'

    trans_table.meta['n_transient'] = n_transient
    trans_table.meta['num_slices'] = 3*num_slices
    trans_table.meta['template'] = transient_template_filename
    trans_table.meta['time_per_slice'] = time_per_slice
    trans_table.meta['bins'] = bins_

    cube_table.write('{}/n{}_s{}_t{}_cube.hdf5'.format(output_path, n_transient, 3*num_slices, transient_template_filename), path='data', overwrite=True)
    trans_table.write('{}/n{}_s{}_t{}_trans.hdf5'.format(output_path, n_transient, 3*num_slices, transient_template_filename), path='data', overwrite=True)


if __name__ == '__main__':
    main()
