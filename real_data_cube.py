import numpy as np
import pandas as pd
import click
import astropy.units as u
from astropy.table import Table
from astropy.time import Time
from fact.io import read_h5py
from IPython import embed


def read_hdf5_file(input_file):
    events = read_h5py(
        input_file,
        key='events',
        columns=[
            'gamma_prediction',
            'ra_prediction',
            'dec_prediction',
            'timestamp'
            ]
    )

    runs = read_h5py(
        input_file,
        key='runs',
        columns=[
            'azimuth',
            'declination',
            'night',
            'ontime',
            'right_ascension',
            'run_start',
            'run_stop',
        ]
    )

    runs = runs.query('ontime > 0').copy()

    return events, runs


def gamma_events_table(events, cut=0.85):
    events = events.query('gamma_prediction > 0.85').copy()
    t = Table(
        {
            'ra': events.ra_prediction.values * u.hourangle,
            'dec': events.dec_prediction.values * u.deg,
            'timestamp': Time(events.timestamp.tolist()),
        }
    )

    return t


def get_runs_per_slice(runs, ontime_max):
    ontime = 0
    slices = [runs.run_start.min()]
    ontime_slice = []
    for i in range(len(runs)):
        run = runs.iloc[i]
    #     print(run)
        ontime += run.ontime
        if ontime > (ontime_max):
            slices.append(run['run_stop'])
            ontime_slice.append(ontime)
            ontime = 0

    return slices, ontime_slice


@click.command()
# @click.argument(
#     'input_file',
#     type=click.Path(file_okay=True, dir_okay=False)
# )
@click.argument(
    'input_file_events',
    type=click.Path(file_okay=True, dir_okay=False)
)
@click.argument(
    'input_file_runs',
    type=click.Path(file_okay=True, dir_okay=False)
)
@click.argument(
    'output_path',
    type=click.Path(file_okay=False, dir_okay=True)
)
@click.option(
    '--source_name',
    '-s',
    type=click.STRING,
    default='Mrk501',
    help='Source Name for output file.'
)
@click.option(
    '--ontime_max',
    '-o',
    type=click.INT,
    default=95*60,
    help='Ontime per slice in seconds. Runs are added to one slice until required ontime is reached.'
)
@click.option(
    '--bins',
    '-b',
    type=click.INT,
    help='Number of ra/dec bins for slice',
    default='80'
)
def main(
    input_file_events,
    input_file_runs,
    # input_file,
    output_path,
    source_name,
    ontime_max,
    bins,
):
    events = pd.read_pickle(input_file_events)
    runs = pd.read_pickle(input_file_runs)

    # events, runs = read_hdf5_file(input_file)

    events_table = gamma_events_table(events)
    slices, ontime_slice = get_runs_per_slice(runs, ontime_max)

    ra_min = events.ra_prediction.min() * u.hourangle.to(u.deg)
    ra_max = events.ra_prediction.max() * u.hourangle.to(u.deg)
    dec_min = events.dec_prediction.min()
    dec_max = events.dec_prediction.max()

    cube = []
    for i in range(len(slices)-1):
        mask = ((events_table['timestamp'] > Time(slices[i])) & (events_table['timestamp'] < Time(slices[i+1])))
        cube.append(
                np.histogram2d(
                    events_table[mask]['ra'].to(u.deg).value,
                    events_table[mask]['dec'].to(u.deg).value,
                    bins=bins
                )[0]
        )
    cube = [cube]
    ontime_slice = [ontime_slice]
    cube_table = Table(
        {
            'cube': cube,
            'ontime': ontime_slice,
        }
    )
    cube_table.meta = {
        'ra': [ra_min, ra_max],
        'dec': [dec_min, dec_max],
        'bins': bins,
        'num_slices': len(slices)
    }
    embed()
    cube_table.write('{}/{}_s{}_cube.hdf5'.format(output_path, source_name, len(slices)), path='data', overwrite=True)


if __name__ == '__main__':
    main()
