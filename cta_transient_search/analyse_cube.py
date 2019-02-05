import pywt
import click
import numpy as np
from tqdm import tqdm
from astropy.table import Table
from ctawave.denoise import thresholding_3d, remove_steady_background


def wavelet_denoising_cube(
    cube_raw,
    n_bg_slices,
    gap,
    bins,
    n_wavelet_slices,
):
    cube_without_steady_source = remove_steady_background(cube_raw, n_bg_slices, gap)
    cube_denoised = [np.zeros([80, 80])]*(gap+n_bg_slices+n_wavelet_slices)

    for i in range(n_wavelet_slices, len(cube_without_steady_source)):

        coeffs = pywt.swtn(
                        data=cube_without_steady_source[i-n_wavelet_slices:i],
                        wavelet='bior1.3',
                        level=2,
                        start_level=0
                    )
        ct = thresholding_3d(coeffs, k=30)
        slice_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')[-1]
        cube_denoised.append(slice_denoised)

    return np.asarray(cube_denoised)


def gradient_benchmark(
    cube_raw
):
    bins = cube_raw.shape[1]
    grad = np.diff(cube_raw, axis=0)
    return np.vstack((np.zeros((1, bins, bins)), grad))


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

    if background:
        cube = cube_raw_table['cube'].data.reshape(bins, bins, -1)
        cube_denoised = gradient_benchmark(cube)
        pos_trigger_pixel = max_pixel_position(cube_denoised)
        list_trigger_position.append(pos_trigger_pixel)
        list_cubes_denoised.append(cube_denoised)

    else:
        for cube in tqdm(cube_raw_table['cube']):
            cube_denoised = gradient_benchmark(cube)
            pos_trigger_pixel = max_pixel_position(cube_denoised)

            list_trigger_position.append(pos_trigger_pixel)
            list_cubes_denoised.append(cube_denoised)

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
        denoised_table.write('{}/n{}_s{}_t{}_bg_gradient_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_cubes,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_bg_gradient_trigger.hdf5'.format(
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
        denoised_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_gradient_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_transient,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                        transient_template_filename,
                                                                        cu_min,
                                                                        z_trans
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_gradient_trigger.hdf5'.format(
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
