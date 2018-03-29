import evaluation
import pandas as pd
import glob
from astropy.table import Table, vstack
import matplotlib.pyplot as plt

from IPython import embed

plt.style.use('ggplot')


def get_metrics_dataframe(path='build/threshold_studies/grid_search/'):
    df_metrics = pd.DataFrame(columns=['Template', 'Threshold', 'TP-Rate', 'FN-Rate', 'FP-Rate', 'TP', 'FN', 'FP', 'num_cubes'])
    for template in range(2, 5):
        input_simulation = Table.read('build/n200_s60_t{}_trans.hdf5'.format(template), path='data')
        for th in range(1, 26):

            input_alert = Table.read('{}n200_s60_t{}_th{}_alert.hdf5'.format(path, template, th), path='data')
            num_cubes = input_alert.meta['n_transient']

            sum_trigger, tp, fp, fn = evaluation.metrics(input_simulation, input_alert)
            df_metrics = df_metrics.append({'Template': template, 'Threshold': th, 'TP-Rate': tp/num_cubes, 'FN-Rate': fn/num_cubes, 'TP': tp, 'FN': fn, 'FP': fp, 'num_cubes': num_cubes}, ignore_index=True)
            # print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))

    return df_metrics


def get_metrics_dataframe_background(path='build/background_studies/grid_search/'):
    df_metrics = pd.DataFrame(columns=['Threshold', 'FP', 'num_cubes'])

    for th in range(1, 26):
        input_alert = Table.read('{}nNone_s60_tNone_th{}_alert.hdf5'.format(path, th), path='data')
        num_cubes = input_alert.meta['n_cubes']

        fp = evaluation.metrics_background(input_alert)
        df_metrics = df_metrics.append({'Threshold': th, 'FP': fp, 'num_cubes': num_cubes}, ignore_index=True)

    return df_metrics


def get_metrics_dataframe_brightness(path='build/threshold_studies/brightness/'):
    df_metrics = pd.DataFrame(columns=['Brightness', 'Threshold', 'TP-Rate', 'FP-Rate', 'TP', 'FP', 'num_cubes'])
    sim_files = glob.glob('{}*_trans.hdf5'.format(path))
    for f in sim_files:
        input_simulation = Table.read(f, path='data')
        for th in range(1, 26):
            input_alert = Table.read('{}n{}_s60_tAll_cu{}_th{}.hdf5'.format(path, input_simulation.meta['n_cubes'], input_simulation.meta['cu_range'], th), path='data')
            num_cubes = input_alert.meta['n_transient']
            sum_trigger, tp, fp, fn = evaluation.metrics(input_simulation, input_alert)
            df_metrics = df_metrics.append({'Brightness': input_simulation.meta['cu_range'], 'Threshold': th, 'TP-Rate': tp/num_cubes, 'FP-Rate': fp/num_cubes, 'TP': tp, 'FN': fn, 'FP': fp, 'num_cubes': tp+fn}, ignore_index=True)
            # print('TP-Rate: {} \n FN-Rate: {} \n FP-Rate: {}'.format(tp/num_cubes, fn/num_cubes, fp/num_cubes))

    return df_metrics


def add_cu_flare(path='build/threshold_studies/grid_search/'):
    for t in range(2, 5):
        trans_table = Table.read('build/n200_s60_t{}_trans.hdf5'.format(t), path='data')
        for th in range(1, 26):
            alert_table = Table.read('{}n200_s60_t{}_th{}_alert.hdf5'.format(path, t, th), path='data')
            alert_table['cu_flare'] = trans_table['cu_flare']
            alert_table.write('{}n200_s60_t{}_th{}_alert.hdf5'.format(path, t, th), path='data', overwrite=True)


def group_alert_tables_by_cu_flare(path='build/threshold_studies/'):
    for i in range(7):
        cu_min = i
        cu_max = i+1
        for th in range(1, 26):
            cu_1_table = Table()
            for t in range(2, 5):
                alert_table = Table.read('{}grid_search/n200_s60_t{}_th{}_alert.hdf5'.format(path, t, th), path='data')
                cu_1_table = vstack([cu_1_table, alert_table[((alert_table['cu_flare'] >= cu_min) & (alert_table['cu_flare'] < cu_max))]])
            cu_1_table.write('{}th_brightness/n{}_s60_tAll_cu{}_th{}.hdf5'.format(path, len(cu_1_table), cu_max, th), path='data', overwrite=True)


def group_trans_tables_by_cu_flare(path='build/threshold_studies/brightness/'):
    for i in range(7):
        cu_min = i
        cu_max = i+1
        cu_1_table = Table()
        for t in range(2, 5):
            trans_table = Table.read('build/n200_s60_t{}_trans.hdf5'.format(t), path='data')
            cu_1_table = vstack([cu_1_table, trans_table[((trans_table['cu_flare'] >= cu_min) & (trans_table['cu_flare'] < cu_max))]])
        cu_1_table.meta['n_cubes'] = len(cu_1_table)
        cu_1_table.meta['cu_range'] = cu_max
        cu_1_table.write('{}th_brightness/n{}_s60_tAll_cu{}_trans.hdf5'.format(path, len(cu_1_table), cu_max), path='data', overwrite=True)


def plot_tp(df, threshold=None):

    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    templates = ['Broad Gaussian', 'Narrow Gaussian', 'Deltapeak + Exponential Decay']
    fig, ax = plt.subplots()
    for t in range(2, 5):
        ax.plot(df[df['Template'] == 2]['Threshold'], df[df['Template'] == t]['TP']/200, label=templates[t-2])
    ax.text(0.95, 0.01, 'PRELIMINARY',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='grey', fontsize=20)
    ax.set_xlabel('Threshold in a. u.')
    ax.set_ylabel('True Positive Rate')
    ax.set_ylim(0, 1)
    ax.legend()
    if threshold is not None:
        ax.axvline(threshold, color='k', linestyle='--')

        plt.savefig('build/plots/true_positive_rate_th_oneW.pdf')
    else:
        plt.savefig('build/plots/true_positive_rate_oneW.pdf')


def plot_fp(df_bg, hours, threshold=None):
    fig, ax = plt.subplots()
    ax.plot(df_bg['Threshold'], df_bg['FP']/hours)
    ax.set_xlabel('Threshold in a. u.')
    ax.set_ylabel('FP per hour')
    ax.text(0.95, 0.04, 'PRELIMINARY',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='grey', fontsize=20)
    if threshold is not None:
        ax.axvline(threshold, color='k', linestyle='--')

        plt.savefig('build/plots/fp_per_hour_th_oneW.pdf')
    else:
        plt.savefig('build/plots/fp_per_hour_oneW.pdf')


def plot_tp_brightness(df):
    fig, ax = plt.subplots()
    for cu in range(1, 8):
        df_cu = df[df['Brightness'] == cu]
        ax.plot(df_cu['Threshold'], df_cu['TP']/df_cu['num_cubes'], label='{}cu'.format(cu))
        ax.set_ylim(0, 1)
        ax.legend(loc=1)
        ax.set_xlabel('Threshold in a. u.')
        ax.set_ylabel('True Positive Rate')
    ax.text(0.95, 0.06, 'PRELIMINARY',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='grey', fontsize=20)
    plt.savefig('build/plots/tp_brightness_oneW.pdf')
    plt.close()


def main():
    # add_cu_flare(path='build/wavelet_studies/grid_search/')
    # group_alert_tables_by_cu_flare(path='build/wavelet_studies/')
    # group_trans_tables_by_cu_flare(path='build/wavelet_studies/')
    df_background = get_metrics_dataframe_background(path='build/wavelet_studies/grid_search/')
    df_signal = get_metrics_dataframe(path='build/wavelet_studies/grid_search/')
    df_brightness = get_metrics_dataframe_brightness(path='build/wavelet_studies/th_brightness/')


    # plot_metrics(df_signal, df_background)
    # plot_accuracy(df_signal, df_background)
    # plot_precision_recall(df_signal, df_background)
    plot_tp_brightness(df_brightness)
    # plot_tp(df_signal, threshold=6)
    plot_tp(df_signal)
    # plot_fp(df_background, 83.33333, threshold=6)
    plot_fp(df_background, 83.33333)



if __name__ == '__main__':
    main()
