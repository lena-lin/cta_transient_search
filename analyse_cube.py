import pywt
import click
import numpy as np
from tqdm import tqdm
from IPython import embed

from astropy.table import Table

from ctawave.denoise import thresholding_3d
from ctawave.toy_models_crab import remove_steady_background


def wavelet_denoising_cube(
    cube_raw,
    n_bg_slices,
    gap,
    bins
):
        if((cube_raw.shape[0] - n_bg_slices - gap) % 4 != 0):
            gap = gap + (cube_raw.shape[0] - n_bg_slices - gap) % 4
        cube = remove_steady_background(cube_raw, n_bg_slices, gap)

        # get wavelet coefficients

        coeffs = pywt.swtn(data=cube, wavelet='bior1.3', level=2, start_level=0)

        # remove noisy coefficents.
        ct = thresholding_3d(coeffs, k=30)
        cube_smoothed = pywt.iswtn(coeffs=ct, wavelet='bior1.3')
        cube_smoothed = np.concatenate([np.zeros([len(cube_raw) - len(cube_smoothed), bins, bins]), cube_smoothed])

        return cube_smoothed


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    'output_path',
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--n_bg_slices',
    '-n_bg',
    help='Number of slices for background mean',
    default=5
)
@click.option(
    '--gap',
    '-g',
    help='Minimal distance to sliding background window (number of slices). Will be enlarged by max. 3 slices if the number of slices for the resulting cube can not be divided by 4.',
    default=5
)
def main(
    input_file,
    output_path,
    n_bg_slices,
    gap,
):
    cube_raw_table = Table.read(input_file, path='data')

    n_transient = cube_raw_table.meta['n_transient']
    num_slices = cube_raw_table.meta['num_slices']
    transient_template_index = cube_raw_table.meta['template']
    bins = cube_raw_table.meta['bins']

    list_cubes_denoised = []
    for cube in tqdm(cube_raw_table['cube']):
        cube_denoised = wavelet_denoising_cube(cube, n_bg_slices, gap, bins)
        list_cubes_denoised.append(cube_denoised)

    denoised_table = Table()
    denoised_table['cube_smoothed'] = list_cubes_denoised
    denoised_table['trans_factor'] = denoised_table['cube_smoothed'].max(axis=2).max(axis=2)

    denoised_table.meta = cube_raw_table.meta
    denoised_table.meta['n_bg_slices'] = n_bg_slices
    denoised_table.meta['gap'] = gap

    denoised_table.write('{}/n{}_s{}_t{}_denoised.hdf5'.format(output_path, n_transient, num_slices, transient_template_index), path='data', overwrite=True)


if __name__ == '__main__':
    main()
