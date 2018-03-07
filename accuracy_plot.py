import evaluation
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt

from IPython import embed

plt.style.use('ggplot')


def get_accuracy_dataframe():
    df_accuracy = pd.DataFrame(columns=['Template', 'Threshold', 'TP-Rate', 'FN-Rate', 'FP-Rate'])
    for template in range(2, 5):
        input_simulation = Table.read('build/n200_s60_t{}_trans.hdf5'.format(template), path='data')
        for th in range(1, 26):
            input_alert = Table.read('build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(template, th), path='data')

            tp, fp, tn, fn, num_cubes = evaluation.accuracy(input_simulation, input_alert)
            df_accuracy = df_accuracy.append({'Template': template, 'Threshold': th, 'TP-Rate': tp/num_cubes, 'FN-Rate': fn/num_cubes, 'FP-Rate': fp/num_cubes}, ignore_index=True)
            # print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))

    return df_accuracy


def get_accuracy_dataframe_background():
    df_accuracy = pd.DataFrame(columns=['Threshold', 'TN-Rate', 'FP-Rate'])

    for th in range(1, 26):
        input_alert = Table.read('build/background_studies/grid_search/nNone_s60_tNone_th{}_alert.hdf5'.format(th), path='data')

        tp, fp, tn, fn, num_cubes = evaluation.accuracy_background(input_alert)
        df_accuracy = df_accuracy.append({'Threshold': th, 'TN-Rate': tn/num_cubes, 'FP-Rate': fp/num_cubes}, ignore_index=True)
        # print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))

    return df_accuracy


def add_cu_flare():
    for t in range(2,5):
        trans_table = Table.read('cta_transient_search/build/n200_s60_t{}_trans.hdf5'.format(t), path='data')
        for th in range(1,26):
            alert_table = Table.read('cta_transient_search/build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(t, th), path='data')
            alert_table['cu_flare'] = trans_table['cu_flare']
            alert_table.write('cta_transient_search/build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(t, th), path='data', overwrite=True)


# def split_accuracy_cu_flare():
#     for t in range(2,5):
#         table =
#         masks = [()]


def plot_accuracy(df):
    fig, ax = plt.subplots()
    ax.plot(df['Threshold'].values, df['FP-Rate'], label='FP-Rate', color='r')
    ax.plot(df['Threshold'].values, df['TN-Rate'], label='TN-Rate', color='g')
    ax.set_ylabel('Rate')
    ax.set_xlabel('Threshold in a.u.')
    ax.set_title('Background')
    ax.legend()
    ax.axvline(4, color='k', linestyle='--')

    plt.savefig('build/plots/accuracy_bg.pdf')
    plt.close()


def main():
    df_accuracy = get_accuracy_dataframe_background()
    plot_accuracy(df_accuracy)


if __name__ == '__main__':
    main()
