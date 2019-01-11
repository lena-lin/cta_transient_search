from astropy.table import Table
import numpy as np

thresholds = [2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]

fp_mean = []
fp_std = []

for th in thresholds:
    table_eval = Table.read('/net/big-tank/POOL/users/llinhoff/CTA/Transients/Alert/nNone_s360000_tNone_th{}_alert.hdf5'.format(th), path='data')

    sum_trigger = []
    for i in range(1000):
        sum_trigger.append(table_eval['trigger_index'][0][i*360 : (i+1) * 360].sum())
    fp_mean.append(np.array(sum_trigger).mean())
    fp_std.append(np.array(sum_trigger).std())
