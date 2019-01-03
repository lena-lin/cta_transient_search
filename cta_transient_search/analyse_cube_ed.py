import pywt
import click
import numpy as np
from tqdm import tqdm
from astropy.table import Table
from ctawave.denoise import thresholding_3d, remove_steady_background
from IPython import embed
import h5py
from collections import deque


def analyze_images(path, n_wavelet_slices, bg_slices, gap):

    size_cube = n_wavelet_slices + bg_slices + gap

    with h5py.File(path) as f:
        data = f['data']
        bins = data.attrs['bins']
        n_images = data.shape[0] // bins // bins

        images = data[:size_cube * bins * bins]['cubes'].reshape(-1, bins, bins)
        queue_bg_sub = deque([])
        cube_denoised = [np.zeros([bins, bins])]*(size_cube - 1)
        for i in range(n_wavelet_slices, n_images):
            queue_bg_sub.append(images[bg_slices + gap - 1] - images[:bg_slices].mean(axis=0))
            if len(queue_bg_sub == n_wavelet_slices):
                coeffs = pywt.swtn(
                    data=np.asarray(queue_bg_sub),
                    wavelet='bior1.3',
                    level=2,
                    start_level=0
                )
                ct = thresholding_3d(coeffs, k=30)
                slice_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')[-1]
                cube_denoised.append(slice_denoised)

                queue_bg_sub.popleft()

        images[:-1] = images[1:]
        images[-1] = data[i * bins * bins: (i+1) * bins * bins]['cubes'].reshape(80, 80)

    return cube_denoised



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
def main(
    input_file,
    output_path,
    n_wavelet_slices,
    n_bg_slices,
    gap,
):
    #cube_raw_table = Table.read(input_file, path='data')
    #slices=np.array(cube_raw_table['cubes']).reshape(360000, 80, 80)
    #cube_raw_table['cubes'] = [slices]

    #n_transient = cube_raw_table.meta.get('n_transient')

    #num_slices = cube_raw_table.meta['num_slices']  # in simulate_cube: 3*n_slices

    #transient_template_index = cube_raw_table.meta.get('template')

    #bins = cube_raw_table.meta['bins']

    list_cubes_denoised = []
    #embed()

    cube_denoised = analyze_images(input_file, 16, 4)
    list_cubes_denoised.append(cube_denoised)

    denoised_table = Table()
    denoised_table['cube_smoothed'] = list_cubes_denoised

    #denoised_table.meta = cube_raw_table.meta
    #denoised_table.meta['n_bg_slices'] = n_bg_slices
    #denoised_table.meta['n_wavelet_slices'] = n_wavelet_slices
    #denoised_table.meta['gap'] = gap

    denoised_table.write('{}/bg_denoised.hdf5'.format(output_path), path='data', overwrite=True)

# Not needed anymore after changes in transient_alert.py
    list_position_c = []

    for c in list_cubes_denoised:
        list_pos = []
        for s in c:
            list_pos.append(np.unravel_index(np.argmax(s, axis=None), s.shape))
        list_position_c.append(list_pos)

    trans_factor_table = Table({'trans_factor': denoised_table['cube_smoothed'].max(axis=2).max(axis=2),
                                'trigger_pos': list_position_c})
    trans_factor_table.meta = denoised_table.meta
    trans_factor_table.write('{}/bg_trigger.hdf5'.format(output_path), path='data', overwrite=True)


if __name__ == '__main__':
    main()
