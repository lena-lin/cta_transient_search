import click
import numpy as np
from astropy.table import Table

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
            distance = np.sqrt(diff_dec**2+deff_ra**2)
            distances.append(distance)
            print(distance)
            if distance <= 1:
                sum_true += 1
            else:
                sum_false += fp
        return sum_true, sum_false , distances

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
    input_simulation, # .trans file
    input_alert, # .alert file
    output_path,
):
        table_simulation = Table.read(input_simulation, path='data')
        table_alert = Table.read(input_alert, path='data')
        num_cubes = table_simulation.meta['n_transient']
        Sum_true, Sum_false, distances = evaluate(table_simulation, table_alert)

        print('True Position: {} \n False Position: {} \n  '.format(Sum_true, Sum_false))
        f = open('{}/evaluation_{}.txt'.format(output_path, num_cubes), 'w')
        f.writelines('Number of simulated transients: {} \n True Position: {} \n False Position: {} \n Distances between predited and true position: {}'.format(num_cubes, Sum_true, Sum_false, distances))
        f.close()


if __name__ == '__main__':
    main()
