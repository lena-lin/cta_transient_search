import click
import numpy as np
from astropy.table import Table

from IPython import embed


def get_next_trigger(trigger_index, start_flare):
    diff_flarestart = [(np.where(trigger_index[i] == True)[0] - start_flare[i]) for i in range(len(trigger_index))]
    closest_trigger = []
    for s in diff_flarestart:
        if len(s) == 0:
            closest_trigger.append(np.nan)
        else:
            closest_trigger.append(s.min())

    return np.asarray(closest_trigger)


def accuracy(table_simulation, table_alert):
    if len(table_simulation) != len(table_alert):
        print('Input tables do not have the same length!')

    else:
        num_cubes = len(table_alert)
        closest_trigger = get_next_trigger(table_alert['trigger_index'], table_simulation['start_flare'])

        tp = np.count_nonzero([abs(closest_trigger) <= 2])
        fn = np.count_nonzero(np.isnan(closest_trigger))
        fp = np.count_nonzero([abs(closest_trigger) > 2])
        tn = 0

        return tp, fp, tn, fn, num_cubes


def accuracy_background(table_alert_bg):
    num_cubes = len(table_alert_bg)
    tn = np.count_nonzero(np.isnan(table_alert_bg['first_trigger'].data))
    fp = np.count_nonzero(~np.isnan(table_alert_bg['first_trigger'].data))
    tp = 0
    fn = 0

    return tp, fp, tn, fn, num_cubes


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
        tp, fp, tn, fn, num_cubes = accuracy(table_simulation, table_alert)

        print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))
        print('TP: {} \n FN: {} \n FP: {} \n TN: {}'.format(tp, fn, fp, tn))
        f = open('build/{}'.format(output_file), 'w')
        f.writelines('Number of simulated transients: {} \n TP: {} \n FN: {} \n FP: {} \n TN: {}'.format(num_cubes, tp, fn, fp, tn))
        f.close()


if __name__ == '__main__':
    main()
