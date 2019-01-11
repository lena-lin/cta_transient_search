from astropy.table import Table
import numpy as np
from evaluation import metrics

thresholds = [2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]

tp_mean = []
tp_std = []

fn_mean = []
fn_std = []

for th in thresholds:
    table_sim = Table.read('/net/big-tank/POOL/users/llinhoff/CTA/Transients/Simulations/n100_s60_t0_cu3.0_trans.hdf5'.format(th), path='data')
    table_alert = Table.read('/net/big-tank/POOL/users/llinhoff/CTA/Transients/Alert/cu3/n100_s60_i0_th{}_alert.hdf5'.format(th), path='data')

    tp = []
    fn = []
    for i in range(10):
        sum_trigger, sum_tp, sum_fp, sum_fn = metrics(table_sim[i*10:(i+1)*10], table_alert[i*10:(i+1)*10])
        tp.append(sum_tp)
        fn.append(sum_fn)

    tp_mean.append(np.array(tp).mean())
    tp_std.append(np.array(tp).std())
    fn_mean.append(np.array(fn).mean())
    fn_std.append(np.array(fn).std())
