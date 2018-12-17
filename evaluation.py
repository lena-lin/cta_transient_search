import click
import numpy as np
from astropy.table import Table

# from IPython import embed


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
    trigger, = np.where(np.diff(ts.astype(int)) == 1)
    if len(trigger) != 0:
        rt = trigger[0]+1    # np.where docu
    # rt = trigger+1  as an array
        print('detected at:', rt, 'simulated at:', start_flare, 'with', np.abs(rt - start_flare))
        tp = np.abs(rt - start_flare) <= 3
        fp = np.abs(rt - start_flare) > 3
        fn = 0
    # if tp == 0:
    #    fn = 1
    # else:
    #    fn = 0
    else:
        print(' not detected ! simulated at:', start_flare)
        fn = 1
        tp = 0
        fp = 0

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


def evaluate(table_simulation, table_alert):
    if len(table_simulation) != len(table_alert):
        print('Input tables do not have the same length! Something went wrong')

    else:
        sum_true = 0
        sum_false = 0
        distances = []
        for prediction, truth in zip(table_alert['pred_position'], table_simulation['position']):
            diff_ra = abs(prediction[0] - truth[0])
            diff_dec = abs(prediction[1] - truth[1])
            distance = np.sqrt(diff_dec**2+diff_ra**2)
            distances.append(distance)
            # print(distance)
            if distance <= 1:  # in degree
                sum_true += 1
            else:
                sum_false += 1
        return sum_true, sum_false, distances


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
    input_simulation,  # Trans
    input_alert,  # Alert
    output_path,
):
        table_simulation = Table.read(input_simulation, path='data')
        Trans_position = table_simulation['position']
        Ra_trans = Trans_position[0][0]  # Each transient at the same Ra/Decs
        Dec_trans = Trans_position[0][1]
        table_alert = Table.read(input_alert, path='data')
        threshold = table_alert.meta['threshold']
        num_cubes = table_simulation.meta['n_transient']
        sum_trigger, tp, fp, fn = metrics(table_simulation, table_alert)

        Sum_true, Sum_false, distances = evaluate(table_simulation, table_alert)

        print('TP: {} \nFN: {} \n FP: {} \n Sum_Trigger: {}'.format(tp, fn, fp, sum_trigger))
        print('Predicted RA in FOV around Crab: ', table_alert['pred_position'][0][0])
        print('Predicted DEC in FOV around Crab: ', table_alert['pred_position'][0][1])
        print('Simulated transient at: ', table_simulation['position'][0][0], table_simulation['position'][0][1])
        print('True Position: {} \n False Position: {} \n  '.format(Sum_true, Sum_false))
        f = open('build/evaluation_{}_{}.txt'.format(num_cubes, threshold), 'w')
        f.writelines('Number of simulated transients: {} \n TP: {} \n FN: {} \n FP: {}'.format(num_cubes, tp, fn, fp))
        f.writelines('\n If Position bool was False, the transient was simulated at RA/DEC =: {}{}'.format(Ra_trans, Dec_trans))
        f.writelines('\n Position evaluation: \n Number true positions: {} \n Number false positions: {} \n Distances between predited and true position for all transients in deg: \n {}'.format(Sum_true, Sum_false, distances))
        f.close()

        return tp, fp, fn


if __name__ == '__main__':
    main()
