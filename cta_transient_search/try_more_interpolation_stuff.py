# coding: utf-8
from astropy.io import fits
irf_path = '/home/lena/Dokumente/CTA/'
irf = fits.open(irf_path)
irf = fits.open('/home/lena/Dokumente/CTA/caldb/data/cta/prod3b/')
irf = fits.open('/home/lena/Dokumente/CTA/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits')
irf.info
irf.info()
irf['BACKGROUND'].header
ipython
irf['BACKGROUND'].data
irf['BACKGROUND'].data['BGD']
plt.imshow(irf['BACKGROUND'].data['BGD'][0])
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib')
plt.imshow(irf['BACKGROUND'].data['BGD'][0])
irf['BACKGROUND'].data['BGD'][0].shape
plt.imshow(irf['BACKGROUND'].data['BGD'][0][0])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-1])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-3])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-5])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-10])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-11])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-8])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-6])
plt.imshow(irf['BACKGROUND'].data['BGD'][0][-5])
irf['BACKGROUND'].data['BGD'][0][-5].min()
irf['BACKGROUND'].data['BGD'][0][-5].max()
bg = irf['BACKGROUND'].data['BGD'][0][-5]
bg[bg>0]
bg[bg>0].min()
bg[bg>0].max() / bg[bg>0].min()
bg = irf['BACKGROUND'].data['BGD'][0][1]
bg[bg>0].max() / bg[bg>0].min()
x = np.arange(bg.shape[0])
import numpy as np
x = np.arange(bg.shape[0])
y = np.arange(bg.shape[1])
x, y = np.meshgrid(x, y)
r = np.sqrt(x**2 + y**2)
plt.plot(r.flatten(), bg.flatten())
get_ipython().magic('matplotlib')
plt.plot(r.flatten(), bg.flatten())
plt.cla()
plt.plot(r.flatten(), bg.flatten(), '.')
x -= len(x) / 2
x = x - (len(x) / 2)
y = y - (len(y) / 2)
r = np.sqrt(x**2 + y**2)
plt.plot(r.flatten(), bg.flatten(), '.')
plt.cla()
plt.plot(r.flatten(), bg.flatten(), '.')
bg = irf['BACKGROUND'].data['BGD'][0][0]
plt.plot(r.flatten(), bg.flatten(), '.')
bg = irf['BACKGROUND'].data['BGD'][0].sum(axis=1)
plt.plot(r.flatten(), bg.flatten(), '.')
bg = irf['BACKGROUND'].data['BGD'][0].sum(axis=0)
plt.plot(r.flatten(), bg.flatten(), '.')
from scipy.interpolate import interp1d
sort_idx = np.argsort(r)
r = r[sort_idx]
bg = bg[sort_idx]
interp1d(bg.cumsum(), r)
len(bg)
interp1d(bg.flatten().cumsum(), r.flatten().cumsum())
f = interp1d(bg.flatten().cumsum(), r.flatten().cumsum())
f(10)
f(1)
f(0)
f(0.0)
f(0.0)
f(0.1)
bg_r = bg.flatten().cumsum() r.flatten().cumsum()
bg_r = bg.flatten().cumsum(); r =  r.flatten()
r = np.append(0, r)
bg_r = np.append(0, bg_r)
f = interp1d(bg_r, r)
f = interp1d(f, bg_r)
f = interp1d(r, bg_r)
f(10)
bg = irf['BACKGROUND'].data['BGD'][0].sum(axis=0)
x = np.arange(bg.shape[0])
x = np.arange(bg.shape[0]) - bg.shape[0] / 2
y = np.arange(bg.shape[1]) - bg.shape[1] / 2
x, y = np.meshgrid(x, y)
x.shape
r = np.sqrt(x **2 + y **2)
bg = bg.flatten()
r.flatten()
plt.plot(r, bg)
r = r.flatten()
plt.plot(r, bg)
plt.cla()
plt.plot(r, bg)
sort_idx = np.argsort(r)
r = r[sort_idx]
bg = bg[sort_idx]
plt.cla()
plt.plot(r, bg)
f = interp1d(r, bg)
test_r = np.linspace(0, 25, 1000)
plt.plot(test_r, f(test_r))
u = np.random.uniform(0, 1, 100000)
f_inv = interp1d(bg.cumsum() / bg.cumsum().max(), r)
random_r = f_inv(u)
r = np.append(0, r)
bg = np.append(0, bg)
f_inv = interp1d(bg.cumsum() / bg.cumsum().max(), r)
random_r = f_inv(u)
plt.hist(random_r, bins=100, range=[0, 25], density=True)
r = _
plt.__version__
import matplotlib
matplotlib.__version__
r[-1].remove()
r[-1].remove()
r[-1]
r
plt.cla()
plt.plot(test_r, f(test_r))
plt.hist(random_r, bins=100, range=[0, 25], normed=True)
f = interp1d(r, bg.cumsum() / bg.cumsum().max())
r
r = np.sqrt(x **2 + y **2)
r = r[sort_idx]
r = r.flatten()[sort_idx]
r = np.append(0, r)
f = interp1d(r, bg.cumsum() / bg.cumsum().max())
F = interp1d(r, bg.cumsum() / bg.cumsum().max())
plt.plot(test_r, F(test_r))
plt.cla()
plt.plot(r, bg / bg.max())
plt.plot(r, bg.cumsum())
plt.plot(r, bg.cumsum() / bg.cumsum())
plt.plot(r, bg.cumsum() / bg.cumsum().max())
plt.cla()
plt.plot(r, bg.cumsum() / bg.cumsum().max())
plt.plot(r, bg / bg.max())
F = interp1d(r, bg.cumsum() / bg.cumsum().max())
plt.plot(test_r, F(test_r))
import pandas as pd
pd.Series(bg)
pd.Series(bg).rolling(5).mean()
plt.plot(r, pd.Series(bg).rolling(5).mean())
plt.plot(r, (pd.Series(bg) / bg.max()).rolling(5).mean())
plt.plot(r, (pd.Series(bg) / bg.max()).rolling(11).mean())
plt.cla()
plt.plot(r, (pd.Series(bg) / bg.max()).rolling(11).mean())
from scipy.interpolate import splrep, splev
spline = splrep(r, bg, xb=0, xe=25, k=3, t=np.linspace(0, 25, 10)[1:-1]) 
spline = splrep(r, bg, k=3, t=np.linspace(0, 25, 10)[1:-1]) 
get_ipython().magic('pinfo splev')
plt.plot(test_r, splev(test_r, spline))
plt.plot(r, bg)
bg[0] = bg[1]
plt.plot(r, bg)
plt.cla()
spline = splrep(r, bg, k=3, t=np.linspace(0, 25, 10)[1:-1]) 
plt.plot(test_r, splev(test_r, spline))
plt.plot(r, bg)
plt.plot(test_r, splev(test_r, spline))
spline = splrep(r, bg.cumsum() / bg.cumsum().max(), k=3, t=np.linspace(0, 25, 10)[1:-1]) 
plt.plot(test_r, splev(test_r, spline))
plt.cla()
plt.plot(test_r, splev(test_r, spline))
y = splev(test_r, spline)
inv_F = interp1d(y, test_r)
u = np.random.uniform(0, 1, 1000)
rand_r = inv_F(u)
t
y
y.min()
r
y
y[0] = 0
r[0]
rand_r = inv_F(u)
inv_F = interp1d(y, test_r)
rand_r = inv_F(u)
plt.hist(rand_r, bins=100, range=[0, 25], normed=True)
plt.plot(r, bg / bg.max())
