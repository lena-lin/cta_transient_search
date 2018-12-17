import matplotlib.pyplot as plt
from matplotlib import animation
from ctawave.plot import TransientPlotter
from astropy.table import Table
import click
plt.style.use('ggplot')


@click.command()
@click.argument('input_cube_raw', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('input_cube_smoothed', type=click.Path(file_okay=True, dir_okay=False))
def main(
    input_cube_raw,
    input_cube_smoothed
):
    cube_raw_table = Table.read(input_cube_raw, path='data')
    cube_denoised_table = Table.read(input_cube_smoothed, path='data')
    cube_with_transient = cube_raw_table['cube'][0]
    cube_smoothed = cube_denoised_table['cube_smoothed'][0]

    time_steps = cube_raw_table.meta['num_slices']
    trans_factor = cube_smoothed.max(axis=1).max(axis=1)
    time_per_slice = cube_raw_table.meta['time_per_slice']

    p = TransientPlotter(cube_with_transient,
                         cube_smoothed,
                         trans_factor,
                         time_per_slice
                         )

    print('Plotting animation. (Be patient)')
    anim = animation.FuncAnimation(
        p.fig,
        p.step,
        #frames=time_steps,
        #interval=5,
        #blit=True,
    )
    from IPython import embed; embed()
    Writer = animation.writers['imagemagick']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('build/gifs/transient_wavelet.gif', writer=writer)


if __name__ == '__main__':
    main()
