import pywt
from scipy.signal import wiener


def remove_steady_background(
            cube_raw,
            n_bg_slices,
            gap
        ):

    if((cube_raw.shape[0] - n_bg_slices - gap) % 4 != 0):
        gap = gap + (cube_raw.shape[0] - n_bg_slices - gap) % 4

    slices = []
    for i in range(n_bg_slices+gap, len(cube_raw)):
        slices.append(cube_raw[i] - cube_raw[(i - gap - n_bg_slices):(i-gap)].mean(axis=0))

    return slices


def remove_steady_background_stationary(
            cube_raw,
            n_bg_slices,
            gap
        ):
    '''
    Remove background by mean over fixed window (without transient)
    to be sure that transient is not substracted from current slice
    '''
    mean_slice = cube_raw[:n_bg_slices+gap].mean(axis=0)

    cube_smoothed = cube_raw - mean_slice

    return cube_smoothed


def remove_steady_source(cube, n_bg_slices, gap=None):
    if gap == None:
        gap = cube.shape[0] - n_bg_slices

    return (cube[-1] - cube[:-(gap + n_bg_slices)].mean(axis=0))


def thresholding_3d(
            coefficient_list,
            sigma_d=2,
            k=3,
            kind='hard',
            sigma_levels=[0.889, 0.2, 0.086, 0.041, 0.020, 0.010, 0.005, 0.0025, 0.0012]):
    '''
    Here we just iterate over all the coefficents and remove those under a certain
    threshold using the pywt.threshold method.
    '''

    r = []
    for level, cs in enumerate(coefficient_list):
        d = {}
        for key, v in cs.items():
            if key == 'aaa':
                d[key] = v
            else:
                d[key] = pywt.threshold(v, sigma_d*k*sigma_levels[level], kind)
        r.append(d)

    return r


def thresholding(
            coefficient_list,
            sigma_d=2,
            k=3,
            kind='hard',
            sigma_levels=[0.889, 0.2, 0.086, 0.041, 0.020, 0.010, 0.005, 0.0025, 0.0012]):
    '''
    at this points we have al coefficients for all planes. We can de-noise by adapting
    small coefficients. Some call this step 'thresholding'.
    When small, or better yet insignificant,
    coefficents are set to 0 without touching the other coefficients
    the process is called hard thresholding.
    Now assume a fixed sigma for the input data noise
    '''
    r = []
    for level, coeffs in enumerate(coefficient_list):
        cA, cs = coeffs
        cs = tuple([pywt.threshold(c, sigma_d*k*sigma_levels[level], kind) for c in cs])
        r.append((cA, cs))
    return r


def wiener_thresholding(coefficient_list):
    '''
    Call the defautl scipy wiener filter on the coefficients for thresholding.
    '''
    r = []
    for level, coeffs in enumerate(coefficient_list):
        cA, cs = coeffs
        cs = tuple([wiener(c) for c in cs])
        r.append((cA, cs))
    return r
