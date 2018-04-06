import pywt
import click
import numpy as np
from tqdm import tqdm
from IPython import embed

from astropy.table import Table

from ctawave.denoise import thresholding_3d, remove_steady_background, remove_steady_background_stationary


def wavelet_denoising_cube(
    cube_raw,
    n_bg_slices,
    gap,
    bins,
    n_wavelet_slices,
):
        cube_without_steady_source = remove_steady_background_stationary(cube_raw, n_bg_slices, gap)

        cube_denoised = []
        for i in range(n_wavelet_slices, len(cube_without_steady_source)):

            coeffs = pywt.swtn(data=cube_without_steady_source[i-n_wavelet_slices:i], wavelet='bior1.3', level=2, start_level=0)
            ct = thresholding_3d(coeffs, k=75)
            slice_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')[-1]
            cube_denoised.append(slice_denoised)

        cube_denoised = np.vstack([np.zeros([n_wavelet_slices, bins, bins]), np.asarray(cube_denoised)])

        return cube_denoised


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--n_wavelet_slices',
    '-w',
    help='Number of slices for wavelet denoising window',
    default=8
)
@click.option(
    '--n_bg_slices',
    '-s',
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
    n_wavelet_slices,
    n_bg_slices,
    gap,
):
    cube_raw_table = Table.read(input_file, path='data')

    try:
        n_transient = cube_raw_table.meta['n_transient']
    except:
        n_transient = None

    num_slices = cube_raw_table.meta['num_slices']

    try:
        transient_template_index = cube_raw_table.meta['template']
    except:
        transient_template_index = None

    bins = cube_raw_table.meta['bins']

    list_cubes_denoised = []
    for cube in tqdm(cube_raw_table['cube']):
        cube_denoised = wavelet_denoising_cube(cube, n_bg_slices, gap, bins, n_wavelet_slices)
        list_cubes_denoised.append(cube_denoised)

    denoised_table = Table()
    denoised_table['cube_smoothed'] = list_cubes_denoised

    denoised_table.meta = cube_raw_table.meta
    denoised_table.meta['n_bg_slices'] = n_bg_slices
    denoised_table.meta['gap'] = gap

    denoised_table.write('{}/n{}_s{}_t{}_w{}_denoised.hdf5'.format(output_path, n_transient, num_slices, transient_template_index, n_wavelet_slices), path='data', overwrite=True)

    trans_factor_table = Table({'trans_factor': denoised_table['cube_smoothed'].max(axis=2).max(axis=2)})
    trans_factor_table.meta = cube_raw_table.meta
    trans_factor_table.write('{}/n{}_s{}_t{}_w{}_trigger.hdf5'.format(output_path, n_transient, num_slices, transient_template_index, n_wavelet_slices), path='data', overwrite=True)


if __name__ == '__main__':
    main()
