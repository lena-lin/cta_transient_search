import click
import numpy as np
from astropy.table import Table

from IPython import embed


def get_next_trigger(trigger_index, start_flare):
    list_trigger = []
    for i in range(len(trigger_index)):
        trigger = trigger_index[i]
        if np.any(trigger):
            list_trigger.append(abs(np.where(trigger)[0] - start_flare[i]).min())
        else:
            list_trigger.append(np.nan)

    return np.asarray(list_trigger)


def count_tp_fp_fn(ts, start_flare):
    rt, = np.where(np.diff(ts.astype(int)) == 1)
    tp = np.any(np.abs(rt - start_flare) <= 3)
    fp = len(rt) - tp
    if tp == 0:
        fn = 1
    else:
        fn = 0

    return tp, fp, fn


def count_fp(ts):
    rt, = np.where(np.diff(ts.astype(int)) == 1)
    fp = len(rt)

    return fp


def metrics(table_simulation, table_alert):
    if len(table_simulation) != len(table_alert):
        print('Input tables do not have the same length!')

    else:
        sum_tp = 0
        sum_fp = 0
        sum_fn = 0
        sum_trigger = 0
        for cube, start_flare in zip(table_alert, table_simulation['start_flare']):
            ts = cube['trigger_index']
            tp, fp, fn = count_tp_fp_fn(ts, start_flare)

            sum_tp += tp
            sum_fp += fp
            sum_fn += fn
            sum_trigger += (tp + fp)

        return sum_trigger, sum_tp, sum_fp, sum_fn


def metrics_background(table_alert_bg):
    sum_fp = 0
    for cube in table_alert_bg:
        ts = cube['trigger_index']
        sum_fp += count_fp(ts)

    return sum_fp



@click.command()
@click.argument('input_simulation', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('input_alert', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
def main(
    input_simulation,
    input_alert,
    output_path,
):
        table_simulation = Table.read(input_simulation, path='data')
        table_alert = Table.read(input_alert, path='data')
        num_cubes = table_simulation.meta['n_transient']
        sum_trigger, tp, fp, fn = metrics(table_simulation, table_alert)

        print('TP: {} \n FN: {} \n FP: {} \n Sum_Trigger: {}'.format(tp, fn, fp, sum_trigger))
        f = open('{}/evaluation_{}.txt'.format(output_path, num_cubes), 'w')
        f.writelines('Number of simulated transients: {} \n TP: {} \n FN: {} \n FP: {}'.format(num_cubes, tp, fn, fp))
        f.close()


if __name__ == '__main__':
    main()
