import pywt
import click
import numpy as np
from tqdm import tqdm
from astropy.table import Table
from ctawave.denoise import thresholding_3d, thresholding, remove_steady_background
from collections import deque


def background_substraction(cube, n_slices_bg, gap_bg):
    if len(cube) != (n_slices_bg + gap_bg + 1):
        print('Invalid cube length of {} for background substraction. Should be {}.'.format(len(cube), (n_slices_bg + gap_bg + 1)))
    else:
        return cube[n_slices_bg + gap_bg] - cube[:n_slices_bg].mean(axis=0)


def background_start(cube, n_slices_bg, gap_bg):
    if len(cube) != (n_slices_bg + gap_bg):
        print('Invalid cube length.')
    else:
        slices_bg = []
        for i in range(1, n_slices_bg + 1):
            slices_bg.append(background_substraction(cube[:i + 1], n_slices_bg=i, gap_bg=0))
        for i in range(1, gap_bg + 1):
            slices_bg.append(background_substraction(cube[:n_slices_bg + i + 1], n_slices_bg, gap_bg=i))

        return slices_bg


def bgSubs_wavelet3d_denoise_lima(cube, n_slices_wavelet, n_slices_off, gap, n_slices_bg, gap_bg):
    alpha = 1 / n_slices_off
    size_bg_cube = n_slices_bg + gap_bg + 1
    queue_bg_sub = deque([])
    queue_denoised = deque([])
    print('Start')
    slices_bg_start = background_start(cube[:n_slices_bg + gap_bg], n_slices_bg, gap_bg)
    print('len start: {}'.format(len(slices_bg_start)))
    for s in slices_bg_start:
        queue_bg_sub.append(s)

    current_slice = n_slices_bg + gap_bg
    while(len(queue_bg_sub) < n_slices_wavelet):
        queue_bg_sub.append(
                            background_substraction(
                                                    cube[(current_slice - size_bg_cube):(current_slice + 1)],
                                                    n_slices_bg,
                                                    gap_bg
                                                    )
                            )
        current_slice += 1

    coeffs_start = pywt.swtn(
                    data=np.array(queue_bg_sub),
                    wavelet='bior1.3',
                    level=2,
                    start_level=0
                )

    ct = thresholding_3d(coeffs_start, k=30)
    start_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')
    for s in start_denoised:
        queue_denoised.append(s)
    print('cube_start finished')

    cube_liMa_S = []
    for i in tqdm(range(current_slice, len(cube))):
        queue_bg_sub.append(
                            background_substraction(
                                                    cube[(current_slice - size_bg_cube):(current_slice + 1)],
                                                    n_slices_bg,
                                                    gap_bg
                                                    )
                            )
        queue_bg_sub.popleft()

        coeffs = pywt.swtn(
                        data=np.array(queue_bg_sub),
                        wavelet='bior1.3',
                        level=2,
                        start_level=0
                    )

        ct = thresholding_3d(coeffs, k=30)
        slice_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')[-1]

        queue_denoised.append(slice_denoised)

        if len(queue_denoised) > (n_slices_off + gap):
            n_off = np.array(queue_denoised)[:n_slices_off].sum(axis=0)
            n_on = slice_denoised
            cube_liMa_S.append(li_ma_significance(n_on, n_off, alpha=alpha))

            queue_denoised.popleft()

    return cube_liMa_S


def wavelet3d_denoise_lima(cube, n_slices_wavelet, n_slices_off, gap):
    alpha = 1 / n_slices_off
    queue_denoised = deque([])
    coeffs_start = pywt.swtn(
                    data=cube[:n_slices_wavelet],
                    wavelet='bior1.3',
                    level=2,
                    start_level=0
                )

    ct = thresholding_3d(coeffs_start, k=30)
    start_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')
    for s in start_denoised:
        queue_denoised.append(s)
    print('cube_start finished')
    cube_liMa_S = []
    for i in tqdm(range(n_slices_wavelet, len(cube))):
        coeffs = pywt.swtn(
                        data=cube[i - n_slices_wavelet + 1:i + 1],
                        wavelet='bior1.3',
                        level=2,
                        start_level=0
                    )

        ct = thresholding_3d(coeffs, k=30)
        slice_denoised = pywt.iswtn(coeffs=ct, wavelet='bior1.3')[-1]

        queue_denoised.append(slice_denoised)

        if len(queue_denoised) > (n_slices_off + gap):
            n_off = np.array(queue_denoised)[:n_slices_off].sum(axis=0)
            n_on = slice_denoised
            cube_liMa_S.append(li_ma_significance(n_on, n_off, alpha=alpha))

            queue_denoised.popleft()

    return cube_liMa_S


def wavelet2d_denoise_lima(cube, n_slices_off, gap):
    alpha = 1 / n_slices_off
    queue_denoised = deque([])
    cube_liMa_S = []
    for slice_raw in cube:
        coeffs = pywt.swt2(
            data=slice_raw,
            wavelet='bior1.3',
            level=2,
            start_level=0
        )
        ct = thresholding(coeffs)
        slice_denoised = pywt.iswt2(coeffs=ct, wavelet='bior1.3')
        queue_denoised.append(slice_denoised)

        if len(queue_denoised) == n_slices_off + gap + 1:
            n_off = np.array(queue_denoised)[:n_slices_off].sum(axis=0)
            n_on = slice_denoised
            cube_liMa_S.append(li_ma_significance(n_on, n_off, alpha=alpha))

            queue_denoised.popleft()

    return cube_liMa_S


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


def li_ma_significance(n_on, n_off, alpha=0.2):
    '''
    Calculate the Li&Ma significance for given
    observations data
    Parameters
    ----------
    n_on: integer or array like
        Number of events for the on observations
    n_off: integer of array like
        Number of events for the off observations
    alpha: float
        Scaling factor for the off observations, for wobble observations
        this is 1 / number of off regions
    '''

    scalar = np.isscalar(n_on)

    n_on = np.array(n_on, copy=False, ndmin=1)
    n_off = np.array(n_off, copy=False, ndmin=1)

    with np.errstate(divide='ignore', invalid='ignore'):
        p_on = n_on / (n_on + n_off)
        p_off = n_off / (n_on + n_off)

        t1 = n_on * np.log(((1 + alpha) / alpha) * p_on)
        t2 = n_off * np.log((1 + alpha) * p_off)

        ts = (t1 + t2)
        significance = np.sqrt(ts * 2)

    significance[np.isnan(significance)] = 0
    significance[n_on < alpha * n_off] = 0

    if scalar:
        return significance[0]

    return significance


def li_ma_benchmark(cube_raw, n_slices_off, gap):
    alpha = 1 / n_slices_off
    slices = []
    for i in range(len(cube_raw) - gap - n_slices_off):
        n_off = cube_raw[i:i + n_slices_off].sum(axis=0)
        n_on = cube_raw[i + n_slices_off + gap]
        slices.append(li_ma_significance(n_on, n_off, alpha=alpha))

    return slices


def gradient_benchmark(
    cube_raw,
    bins
):
    grad = np.diff(cube_raw, axis=0)
    # return np.vstack((np.zeros((1, bins, bins)), grad))
    return grad


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
        cube_S = bgSubs_wavelet3d_denoise_lima(cube, 8, 5, 3, 3, 3)
        pos_trigger_pixel = max_pixel_position(cube_S)
        list_trigger_position.append(pos_trigger_pixel)
        list_cubes_denoised.append(cube_S)

    else:
        print('signal')
        for cube in tqdm(cube_raw_table['cube']):
            cube_S = bgSubs_wavelet3d_denoise_lima(cube, 8, 5, 3, 3, 3)
            pos_trigger_pixel = max_pixel_position(cube_S)

            list_trigger_position.append(pos_trigger_pixel)
            list_cubes_denoised.append(cube_S)

    from IPython import embed; embed()

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
        denoised_table.write('{}/n{}_s{}_t{}_bg_bgSubs_3dwavelet_lima_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_cubes,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_bg_bgSubs_3dwavelet_lima_trigger.hdf5'.format(
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
        denoised_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_bgSubs_3dwavelet_lima_denoised.hdf5'.format(
                                                                        output_path,
                                                                        n_transient,
                                                                        num_slices,
                                                                        time_per_slice,
                                                                        transient_template_filename,
                                                                        cu_min,
                                                                        z_trans
                                                                    ), path='data', overwrite=True)

        trans_factor_table.meta = denoised_table.meta
        trans_factor_table.write('{}/n{}_s{}_t{}_i{}_cu{}_z{}_bgSubs_3dwavelet_lima_trigger.hdf5'.format(
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
