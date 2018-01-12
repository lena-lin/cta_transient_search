import click
import numpy as np
from astropy.table import Table

from IPython import embed

@click.command()
@click.argument('input_simulation', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('input_alert', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('output_file', type=click.Path(file_okay=True, dir_okay=False))

def main(
    input_simulation,
    input_alert,
    output_file
    ):


    table_simulation = Table.read(input_simulation, path='data')
    table_alert = Table.read(input_alert, path='data')

    if len(table_simulation) != len(table_alert):
        print('Input tables do not have the same length!')

    else:
        num_cubes = len(table_alert)
        diff_flarestart = np.asarray(table_simulation['start_flare'] - table_alert['first_trigger'])
        tp = len(diff_flarestart[diff_flarestart == 0])
        fn = len(diff_flarestart[np.isnan(diff_flarestart)])
        fp = len(diff_flarestart[diff_flarestart != 0]) - fn
        tn = 0

        embed()

        print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))


if __name__ == '__main__':
    main()
