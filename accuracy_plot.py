import evaluation
import pandas as pd
import glob
from astropy.table import Table, vstack
import matplotlib.pyplot as plt

from IPython import embed

plt.style.use('ggplot')


def get_accuracy_dataframe():
    df_accuracy = pd.DataFrame(columns=['Template', 'Threshold', 'TP-Rate', 'FN-Rate', 'FP-Rate', 'num_cubes'])
    for template in range(2, 5):
        input_simulation = Table.read('build/n200_s60_t{}_trans.hdf5'.format(template), path='data')
        for th in range(1, 26):
            input_alert = Table.read('build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(template, th), path='data')

            tp, fp, tn, fn, num_cubes = evaluation.accuracy(input_simulation, input_alert)
            df_accuracy = df_accuracy.append({'Template': template, 'Threshold': th, 'TP-Rate': tp/num_cubes, 'FN-Rate': fn/num_cubes, 'FP-Rate': fp/num_cubes, 'num_cubes': num_cubes}, ignore_index=True)
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


def get_accuracy_dataframe_brightness():
    df_accuracy = pd.DataFrame(columns=['Brightness', 'Threshold', 'TP-Rate', 'FP-Rate', 'TP', 'FP', 'num_cubes'])
    sim_files = glob.glob('build/threshold_studies/brightness/*_trans.hdf5')
    for f in sim_files:
        input_simulation = Table.read(f, path='data')
        for th in range(1, 26):
            input_alert = Table.read('build/threshold_studies/brightness/n{}_s60_tAll_cu{}_th{}.hdf5'.format(input_simulation.meta['n_cubes'], input_simulation.meta['cu_range'], th), path='data')

            tp, fp, tn, fn, num_cubes = evaluation.accuracy(input_simulation, input_alert)
            df_accuracy = df_accuracy.append({'Brightness': input_simulation.meta['cu_range'], 'Threshold': th, 'TP-Rate': tp/num_cubes, 'FP-Rate': fp/num_cubes, 'TP': tp, 'FP': fn, 'num_cubes': num_cubes}, ignore_index=True)
            # print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))

    return df_accuracy


def add_cu_flare():
    for t in range(2, 5):
        trans_table = Table.read('cta_transient_search/build/n200_s60_t{}_trans.hdf5'.format(t), path='data')
        for th in range(1, 26):
            alert_table = Table.read('cta_transient_search/build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(t, th), path='data')
            alert_table['cu_flare'] = trans_table['cu_flare']
            alert_table.write('cta_transient_search/build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(t, th), path='data', overwrite=True)


def group_alert_tables_by_cu_flare():
    for i in range(7):
        cu_min = i
        cu_max = i+1
        for th in range(1,26):
            cu_1_table = Table()
            for t in range(3,5):
                alert_table = Table.read('cta_transient_search/build/threshold_studies/grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(t, th), path='data')
                cu_1_table = vstack([cu_1_table, alert_table[((alert_table['cu_flare'] >= cu_min) & (alert_table['cu_flare'] < cu_max))]])
            cu_1_table.write('cta_transient_search/build/threshold_studies/brightness/n{}_s60_tAll_cu{}_th{}.hdf5'.format(len(cu_1_table), cu_max, th), path='data', overwrite=True)


def group_trans_tables_by_cu_flare():
    for i in range(7):
        cu_min = i
        cu_max = i+1
        cu_1_table = Table()
        for t in range(2, 5):
            trans_table = Table.read('cta_transient_search/build/n200_s60_t{}_trans.hdf5'.format(t), path='data')
            cu_1_table = vstack([cu_1_table, trans_table[((trans_table['cu_flare'] >= cu_min) & (trans_table['cu_flare'] < cu_max))]])
        cu_1_table.meta['n_cubes'] = len(cu_1_table)
        cu_1_table.meta['cu_range'] = cu_max
        cu_1_table.write('cta_transient_search/build/threshold_studies/brightness/n{}_s60_tAll_cu{}_trans.hdf5'.format(len(cu_1_table), cu_max), path='data', overwrite=True)

# def split_accuracy_cu_flare():
#     for t in range(2,5):
#         table =
#         masks = [()]


def plot_accuracy_brightness(df_fn_cu, df_fp):
    fig, ax = plt.subplots()
    for cu in range(1, 8):
        ax.plot(df[df['Brightness'] == cu]['FP-Rate'], df[df['Brightness'] == cu]['TP-Rate'], label='{}-{}'.format(cu-1, cu))
    ax.set_ylabel('TP-Rate')
    ax.set_xlabel('FP-Rate')
    ax.set_title('ROC')
    ax.legend()

    plt.savefig('build/plots/accuracy_brightness.pdf')
    plt.close()


def plot_accuracy_background(df):
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


def plot_accuracy(df, df_fp):
    templates = ['Broad Gaussian', 'Narrow Gaussian', 'Deltapeak + Exponential Decay']
    for t in range(2, 5):
        fig, ax = plt.subplots()
        #ax.plot(df[df['Template'] == t]['Threshold'].values, df[df['Template'] == t]['FP-Rate'], label='FP-Rate', color='r')
        ax.plot(df_fp['Threshold'].values, df_fp['FP-Rate'], label='FP-Rate', color='r')
        ax.plot(df[df['Template'] == t]['Threshold'].values, df[df['Template'] == t]['TP-Rate'], label='TP-Rate', color='g')
        ax.plot(df[df['Template'] == t]['Threshold'].values, df[df['Template'] == t]['FN-Rate'], label='FN-Rate', color='b')
        ax.set_ylabel('Rate')
        ax.set_xlabel('Threshold in a.u.')
        ax.set_title(templates[t-2])
        ax.legend()
        ax.axvline(4, color='k', linestyle='--')

        plt.savefig('build/plots/accuracy_t{}.pdf'.format(templates[t-2]))
        plt.close()


def main():
    df_accuracy_fp = get_accuracy_dataframe_background()
    df_accuracy = get_accuracy_dataframe()

    plot_accuracy(df_accuracy, df_accuracy_fp)


if __name__ == '__main__':
    main()
