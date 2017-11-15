import matplotlib.pyplot as plt
import pywt
import click
import numpy as np


from astropy.io import fits

from ctawave.denoise import thresholding_3d
from ctawave.toy_models_crab import remove_steady_background
plt.style.use('ggplot')


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('out_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '--time_per_slice',
    '-t',
    type=click.INT,
    help='Measuring time for one slice in s',
    default='30'
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
@click.option(
    '--cmap',
    '-c',
    help='Colormap to use for histograms',
    default='viridis'
)
@click.option(
    '--bins',
    '-b',
    help='Colormap to use for histograms',
    default=80
)


def main(
    input_file,
    out_file,
    time_per_slice,
    n_bg_slices,
    gap,
    cmap,
    bins
        ):
        cube_raw = fits.open(input_file)[0].data.reshape([-1, bins, bins])
        if((cube_raw.shape[0] - n_bg_slices - gap) % 4 != 0):
            gap = gap + (cube_raw.shape[0] - n_bg_slices - gap) % 4
        cube = remove_steady_background(cube_raw, n_bg_slices, gap, [bins, bins])

        # get wavelet coefficients

        coeffs = pywt.swtn(data=cube, wavelet='bior1.3', level=2, start_level=0)

        # remove noisy coefficents.
        ct = thresholding_3d(coeffs, k=30)
        cube_smoothed = pywt.iswtn(coeffs=ct, wavelet='bior1.3')
        # embed()
        cube_smoothed = np.concatenate([np.zeros([len(cube_raw) - len(cube_smoothed), bins, bins]), cube_smoothed])

        # some Criterion which could be used to trigger this.
        trans_factor = cube_smoothed.max(axis=1).max(axis=1)

        # save trans_factor values for further evaluation (fits)
        hdu_trans_factor = fits.PrimaryHDU(trans_factor)
        hdu_trans_factor.writeto(out_file, overwrite=True)


if __name__ == '__main__':
    main()
