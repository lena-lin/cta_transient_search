import click
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

# from IPython import embed

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
    table = input_file
    table['trans_factor_mov_avg'] = list(map(moving_average, table['trans_factor']))
    table['trans_factor_diff'] = table['trans_factor'] - table['trans_factor_mov_avg']

    return table


def get_next_trigger(trigger_index, start_flare):
    list_trigger = []
    for i in range(len(trigger_index)):
        trigger = trigger_index[i]
        if np.any(trigger):
            list_trigger.append(abs(np.where(trigger)[0] - start_flare[i]).min())
        else:
            list_trigger.append(np.nan)

    return np.asarray(list_trigger)


def send_alert(table, threshold):
    trigger_index = table['trans_factor_diff'] > threshold
    found_trigger = trigger_index.sum(axis=1)

    return trigger_index, found_trigger


def get_transient_position(
    list_cubes,
    first_trigger,
    fov,
    bins,
    source,
):
    #source_coordinates = SkyCoord.from_name(source)
    source_coordinates = SkyCoord('05 34 31.97 +22 00 52.1', unit=(u.hourangle, u.deg))
    list_positions = []
    for trig, cube in zip(first_trigger, list_cubes):
        if len(trig[trig!=False]) > 0:
            trigger = np.where(np.diff(trig==1))[0][0]+1  # simple test with np.where
            slice = cube[trigger]
            max_pos = np.unravel_index(np.argmax(slice), slice.shape)
            max_pos_ra = max_pos[0] * fov/bins + source_coordinates.ra.deg - fov/2
            max_pos_dec = max_pos[1] * fov/bins + source_coordinates.dec.deg - fov/2
            #max_pos_ra = np.interp(max_pos[0],[0,bins],[source_coordinates.ra.deg - fov / 2,source_coordinates.ra.deg + fov / 2])
            #max_pos_ra = np.interp(max_pos[1],[0,bins],[source_coordinates.dec.deg - fov / 2,source_coordinates.dec.deg + fov / 2])
        else:
            max_pos_ra = np.nan
            max_pos_dec = np.nan
        list_positions.append([max_pos_ra, max_pos_dec])

    return list_positions


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '--output_path',
    type=click.Path(dir_okay=True),
    help='Directory for output file (astropy table)',
    default='build'
)
@click.option(
    '--threshold',
    '-t',
    help='Trigger threshold, default=6',
    default=4
)
def main(
    input_file,  # ! .trigger nicht denoised table
    output_path,
    threshold
):
    table_den = Table.read(input_file, path='data')
    trans_factor_table = Table({'trans_factor': table_den['cube_smoothed'].max(axis=2).max(axis=2)})
    trans_factor_table.meta = table_den.meta
    denoised_table = get_smoothed_table(trans_factor_table)
    trigger_index, found_trigger = send_alert(denoised_table, threshold)


    try:
        n_transient = denoised_table.meta['n_transient']
    except:
        n_transient = None

    num_slices = denoised_table.meta['num_slices']

    try:
        transient_template_index = denoised_table.meta['template']
    except:
        transient_template_index = None

    num_slices = denoised_table.meta['num_slices']

    alert_table = Table()
    alert_table['trigger_index'] = trigger_index  # list of bools (len=number slices), true for trigger, false for no trigger
    alert_table['found_trigger'] = found_trigger  # number of triggers found in series (aka number of true in trigger index)
    alert_table['trans_factor_diff'] = denoised_table['trans_factor_diff']  # time trigger criterion

    alert_table['pred_position'] = get_transient_position(
                                         table_den['cube_smoothed'],
                                         trigger_index, denoised_table.meta['fov'],
                                         denoised_table.meta['bins'],
                                         denoised_table.meta['steady_source']
                                     )

    alert_table.meta = denoised_table.meta
    alert_table.meta['threshold'] = threshold
    #alert_table.write('{}/n{}_s{}_t{}_th{}_alert.hdf5'.format(output_path, n_transient, num_slices, transient_template_index, threshold), path='data', overwrite=True)
    alert_table.write('{}/n{}_s{}_t{}_alert.hdf5'.format(output_path, n_transient, num_slices, transient_template_index), path='data', overwrite=True)


if __name__ == '__main__':
    main()
