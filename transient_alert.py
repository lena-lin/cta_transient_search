import click
import numpy as np
from astropy.table import Table

'''
    Send alerts by analysing the timeseries from analyse_cube.py
    Input: astropy table (hdf5) containing one timeseries per simulation of length num_slices (60)
    Output: astropy table with 3 columns: trigger_index (list of num_slices bools, True: trans_factor over threshold),
            found_trigger (number of values above threshold (trigger)), first_trigger (index of first trigger slice)
'''


def moving_average(timeseries, interval_size=10):
    list_averages = np.zeros(interval_size).tolist()
    for i in range(interval_size, len(timeseries)):
        list_averages.append(timeseries[i-interval_size:i].mean())
    averages = np.asarray(list_averages)

    return averages


def get_smoothed_table(input_file):
    table = Table.read(input_file, path='data')
    table['trans_factor_mov_avg'] = list(map(moving_average, table['trans_factor']))
    table['trans_factor_diff'] = table['trans_factor'] - table['trans_factor_mov_avg']

    return table


def get_th_index(timeseries, threshold):
    return [i for i, diff in enumerate(timeseries) if diff > threshold]


def send_alert(table, threshold):
    trigger_index = table['trans_factor_diff'] > threshold
    found_trigger = trigger_index.sum(axis=1)
    list_trigger = [np.where(series == True)[0] for series in trigger_index]

    first_trigger = []
    for e in list_trigger:
        if len(e) > 0:
            first_trigger.append(e[0])
        else:
            first_trigger.append(np.nan)

    return trigger_index, found_trigger, first_trigger


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    'output_path',
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--threshold',
    '-th',
    default=6
)
def main(
    input_file,
    output_path,
    threshold
):
    timeseries_table = get_smoothed_table(input_file)
    trigger_index, found_trigger, first_trigger = send_alert(timeseries_table, threshold)

    n_transient = timeseries_table.meta['n_transient']
    num_slices = timeseries_table.meta['num_slices']
    transient_template_index = timeseries_table.meta['template']

    alert_table = Table()
    alert_table['trigger_index'] = trigger_index  # list of bools (len=number slices), true for trigger, false for no trigger
    alert_table['found_trigger'] = found_trigger  # number of triggers found in series (aka number of true in trigger index)
    alert_table['first_trigger'] = first_trigger  # first trigger slice

    alert_table.meta = timeseries_table.meta
    alert_table.meta['threshold'] = threshold

    alert_table.write('{}/n{}_s{}_t{}_alert.hdf5'.format(output_path, n_transient, num_slices, transient_template_index), path='data', overwrite=True)


if __name__ == '__main__':
    main()
