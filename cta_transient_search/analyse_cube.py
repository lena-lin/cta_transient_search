import pywt
import click
import numpy as np
from tqdm import tqdm
from astropy.table import Table
from ctawave.denoise import thresholding_3d, remove_steady_background
from scipy import special
from skimage.filters import gaussian
from ctawave.denoise import thresholding_3d, thresholding, remove_steady_background
from collections import deque


def smooth_psf_kernel(cube_raw, psf=0.1):
    bins = cube[0].shape[0]
    fov = 8
    list_denoised = []
    for slice in cube_raw:
        denoised_slice = gaussian(slice, sigma=psf / (fov / bins))
        list_denoised.append(denoised_slice)
    return list_denoised


def denoise_slice(cube, n_slices_bg, gap_bg, psf=0.1):
    if len(cube) != (n_slices_bg + gap_bg + 1):
        print('Invalid cube length of {} for background substraction. Should be {}.'.format(len(cube), (n_slices_bg + gap_bg + 1)))
    else:
        bins = cube[0].shape[0]
        fov = 8
        bg_subs_slice = cube[n_slices_bg + gap_bg] - cube[:n_slices_bg].mean(axis=0)
        smoothed_slice = gaussian(bg_subs_slice, sigma=psf / (fov / bins))

        return smoothed_slice


def denoise_start(cube, n_slices_bg, gap_bg):
    if len(cube) != (n_slices_bg + gap_bg):
        print('Invalid cube length.')
    else:
        slices_denoised = []
        for i in range(1, n_slices_bg + 1):
            # print('n_slices_bg: {}'.format(i))
            slices_denoised.append(denoise_slice(cube[:i + 1], n_slices_bg=i, gap_bg=0))
        for k in range(1, gap_bg):
            # print('n_slices_bg: {}'.format(k))
            slices_denoised.append(denoise_slice(cube[:n_slices_bg + k + 1], n_slices_bg, gap_bg=k))

        return slices_denoised


def bgSubs_gaussian_smooth(cube, n_slices_bg, gap_bg):
    size_bg_cube = n_slices_bg + gap_bg + 1
    slices_denoised = denoise_start(cube[:n_slices_bg + gap_bg], n_slices_bg, gap_bg)

    for i in range(size_bg_cube, len(cube)):
        slices_denoised.append(denoise_slice(cube[(i - size_bg_cube + 1):(i + 1)], n_slices_bg, gap_bg))

    return slices_denoised


def max_pixel_position(
    cube_denoised,
):
    list_pos = []
    for s in cube_denoised:
        list_pos.append(np.unravel_index(np.argmax(s, axis=None), s.shape))
    return list_pos


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
    default=4
)
@click.option(
    '--gap',
    '-g',
    help='Minimal distance to sliding background window (number of slices). Will be enlarged by max. 3 slices if the number of slices for the resulting cube can not be divided by 4.',
    default=4
)
@click.option(
    '--background',
    '-b',
    help='Boolean set True if background data without transient signal is analyzed',
    is_flag=True
)
def main(
    input_file,
    output_path,
    n_wavelet_slices,
    n_bg_slices,
    gap,
    background,
):
    cube_raw_table = Table.read(input_file, path='data')
    bins = cube_raw_table.meta['bins']

    list_cubes_denoised = []
    list_trigger_position = []

    if background == True:
        cube = cube_raw_table['cube'].data.reshape(-1, bins, bins)
        print('bg', cube.shape)
        cube_smoothed = bgSubs_gaussian_smooth(cube, 5, 3)
        pos_trigger_pixel = max_pixel_position(cube_smoothed)
        list_trigger_position.append(pos_trigger_pixel)
        list_cubes_denoised.append(cube_S)

    else:
        print('signal')
        for cube in tqdm(cube_raw_table['cube']):
            cube_smoothed = bgSubs_gaussian_smooth(cube, 5, 3)
            pos_trigger_pixel = max_pixel_position(cube_smoothed)

            list_trigger_position.append(pos_trigger_pixel)
            list_cubes_denoised.append(cube_smoothed)

    denoised_table = Table()
    denoised_table['cube_smoothed'] = list_cubes_denoised

    denoised_table.meta = cube_raw_table.meta
    denoised_table.meta['n_bg_slices'] = n_bg_slices
    denoised_table.meta['n_wavelet_slices'] = n_wavelet_slices
    denoised_table.meta['gap'] = gap

    trans_factor_table = Table({'trans_factor': denoised_table['cube_smoothed'].max(axis=2).max(axis=2),
                                'trigger_pos': list_trigger_position})

    if background:
        num_slices = cube_raw_table.meta['num_slices']  # in simulate_cube: 3*n_slices
        time_per_slice = cube_raw_table.meta['time_per_slice']
        n_cubes = cube_raw_table.meta['n_cubes']
        denoised_table.write('{}/n{}_s{}_t{}_bg_bgsubs_smooth_gaussian_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_cubes,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_bg_bgsubs_smooth_gaussian_trigger.hdf5'.format(
                                                                        output_path,
                                                                        n_cubes,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                    ), path='data', overwrite=True)
    else:
        num_slices = cube_raw_table.meta['num_slices']  # in simulate_cube: 3*n_slices
        time_per_slice = cube_raw_table.meta['time_per_slice']
        n_transient = cube_raw_table.meta['n_transient']
        transient_template_filename = cube_raw_table.meta['template']
        cu_min = cube_raw_table.meta['min brightness in cu']
        z_trans = cube_raw_table.meta['redshift']
        denoised_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_bgsubs_smooth_gaussian_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_transient,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                        transient_template_filename,
                                                                        cu_min,
                                                                        z_trans
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_bgsubs_smooth_gaussian_trigger.hdf5'.format(
                                                                        output_path,
                                                                        n_transient,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                        transient_template_filename,
                                                                        cu_min,
                                                                        z_trans
                                                                    ), path='data', overwrite=True)


if __name__ == '__main__':
    main()
