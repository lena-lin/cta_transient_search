import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
plt.style.use('ggplot')
import numpy as np

class TransientPlotter(object):

    @staticmethod
    def plot_trigger_criterion(transient):
        fig, ax = plt.subplots(1)
        # rotate and align the tick labels so they look better
        for i in range(len(transient.trigger_criterion)):
            ax.plot(transient.trigger_criterion_timestamps[i], transient.trigger_criterion[i], '.')

        import matplotlib.dates as mdates
        ax.fmt_xdata = mdates.DateFormatter('%H:%M:%S')
        fig.autofmt_xdate()

        return fig

    def __init__(self, left_cube, right_cube, trans_factor, trans_scale, time_per_slice, cmap='viridis'):
        self.fig = plt.figure(figsize=(12, 10))

        ax1 = plt.subplot2grid((2, 2), (0, 0))
        ax2 = plt.subplot2grid((2, 2), (0, 1))
        ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
        # ax4 = plt.subplot2grid((3, 2), (1, 0), colspan=2)

        ax1.tick_params(labelbottom='off', labelleft='off')
        ax2.tick_params(labelbottom='off', labelleft='off')

        # ax3.set_xticks(ax3.get_xticks() * time_per_slice)
        ax3.set_xlabel('Time Step in a.u.')
        ax3.set_ylabel('Trigger Criterion in a.u.')

        # ax4.set_ylabel('Transient Scale Factor in a.u.')

        vmax = left_cube.max()
        self.l_quad = ax1.pcolormesh(left_cube[0], cmap=cmap, vmin=0, vmax=vmax/5)
        self.r_quad = ax2.pcolormesh(left_cube[0], cmap=cmap, vmin=0, vmax=vmax/5)

        self.line,  = ax3.plot(0, trans_factor[0])
        # self.trans,  = ax4.plot(0, trans_scale[0])

        ax3.set_xlim([0, len(trans_factor)])
        ax3.set_ylim([0, trans_factor.max() + 1])

        # ax4.set_xlim([0, len(trans_scale)])
        # ax4.set_ylim([0, trans_scale.max() + 0.1])

        self.left_cube = left_cube
        self.right_cube = right_cube
        self.trans_factor = trans_factor
        # self.trans_scale = trans_scale
        self.x = []
        self.y = []
        # self.y4 = []

    def step(self, t):
        self.x.append(t)
        self.y.append(self.trans_factor[t])
        # self.y4.append(self.trans_scale[t])

        l = self.left_cube[t]
        r = self.right_cube[t]
        self.l_quad.set_array(l.ravel())
        self.r_quad.set_array(r.ravel())
        self.line.set_data(self.x, self.y)
        # self.trans.set_data(self.x, self.y4)

        return [self.l_quad, self.r_quad, self.line]


class CubePlotter(object):

    def __init__(self, cube, source, fov=12 * u.deg, vmax=None, cmap='viridis'):
        crab_coord = SkyCoord.from_name(source)
        fig, ax = plt.subplots(1, 1)
        self.fig = fig

        ra_ticks = np.linspace(crab_coord.ra.deg - fov.value / 2, crab_coord.ra.deg + fov.value / 2, cube.shape[0])
        dec_ticks = np.linspace(crab_coord.dec.deg - fov.value / 2, crab_coord.dec.deg + fov.value / 2, cube.shape[1])
        X, Y = np.meshgrid(ra_ticks, dec_ticks)

        if vmax == None:
            vmax = cube.max()

        ax.set_aspect(1)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        self.quad = ax.pcolormesh(X, Y, cube, cmap=cmap, vmin=0, vmax=vmax)

        self.cube = cube

    def step(self, t):
        l = self.cube[t]

        self.quad.set_array(l.ravel())
        return [self.quad]




def pixel_histogram(*images, labels=['reconstructed pixel values'], bins=40):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for image, label in zip(images, labels):
        ax.hist(image.flatten(), bins=bins, histtype='step', label=label)

    plt.legend()
    return fig

def coefficients(coeff_list, cmap='gray'):
    titles = ['Approximation', ' Horizontal detail',
              'Vertical detail', 'Diagonal detail']
    fig, axs = plt.subplots(nrows=len(coeff_list), ncols=4)
    for i, coeffs in enumerate(coeff_list):
        cA, (cH, cV, cD) = coeffs
        axes_row = axs[i]
        for ax, c in zip(axes_row, [cA, cH, cV, cD]):
            im = ax.imshow(c, origin='image', interpolation="nearest", cmap=cmap)
            ax.tick_params(labelbottom='off', labelleft='off')
            cbar = fig.colorbar(im, ax=ax)
            cbar.ax.tick_params(labelsize=6)
            # ax.tick_params(axis='both', which='both', labelsize=4)

    fig.suptitle("swt2 coefficients", fontsize=12)
    fig.tight_layout()

    return fig

def results(original_image, reconstructed_image, npe_truth, cmap='viridis'):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    im = ax1.imshow(original_image, interpolation='nearest', cmap=cmap)
    fig.colorbar(im, ax=ax1)
    ax1.set_title('original')
    ax1.tick_params(labelbottom='off', labelleft='off')

    im = ax2.imshow(reconstructed_image, interpolation='nearest', cmap=cmap)
    fig.colorbar(im, ax=ax2)
    ax2.set_title('reconstructed')
    ax2.tick_params(labelbottom='off', labelleft='off')

    im = ax3.imshow(original_image - reconstructed_image, interpolation='nearest', cmap=cmap)
    fig.colorbar(im, ax=ax3)
    ax3.set_title('residual')
    ax3.tick_params(labelbottom='off', labelleft='off')

    im = ax4.imshow(npe_truth, interpolation='nearest', cmap=cmap)
    fig.colorbar(im, ax=ax4)
    ax4.set_title('simulated photons')
    ax4.tick_params(labelbottom='off', labelleft='off')

    return fig
