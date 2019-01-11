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
%matplotlib
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
%matplotlib
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
splev?
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
%save try_more_interpolation_stuff.py 1-180
spline
r
bg
bg
len(bg)
len(r)
36*36
len(np.unique(r))
df_bg = pd.DataFrame(data={'r': r, 'bg': bg})
df_bg
df_bg.groupby(['r']).mean()
df_bg_grouped = df_bg.groupby(['r']).mean()
plt.plot(df_bg_grouped['r'], df_bg_grouped['bg'])
df_bg_grouped
df_bg_grouped.index
df_bg_grouped.index[0]
df_bg_grouped.index.values
df_bg_grouped['bg']
df_bg_grouped['bg'].values
plt.plot(df_bg_grouped.index.values, df_bg_grouped['bg'].values)
plt.cla()
plt.plot(df_bg_grouped.index.values, df_bg_grouped['bg'].values)
plt.plot(df_bg_grouped.index.values, df_bg_grouped['bg'].values, 'o')
plt.plot(r, bg)
r_unique = df_bg_grouped.index.values; bg_mean =  df_bg_grouped['bg'].values)
r_unique = df_bg_grouped.index.values; bg_mean =  df_bg_grouped['bg'].values
spline = splrep(r_unique, bg_mean.cumsum() / bg_mean.cumsum().max(), k=3, t=np.linspace(0, 25, 10)[1:-1])
spline
plt.plot(r_unique, bg_mean.cumsum() / bg_mean.max())
plt.plot(r_unique, bg_mean.cumsum() / bg_mean.cumsum().max())
plt.cla()
plt.plot(r_unique, bg_mean.cumsum() / bg_mean.cumsum().max())
plt.plot(r, bg)
plt.plot(df_bg_grouped.index.values, df_bg_grouped['bg'].values, 'o')
spline
splev(r)
splev(test_r, spline)
plt.plot(test_r, splev(test_r, spline))
spline = splrep(bg_mean.cumsum() / bg_mean.cumsum().max(), r_unique, k=3, t=np.linspace(0, 1, 10)[1:-1])
plt.plot(test_r, splev(np.random.randint(0,1,50), spline))
np.random.randint(0,1,50)
np.random.rand(0,1,50)
np.random.uniform(50)
np.random.uniform(0,1,50)
test_u = np.random.uniform(0,1,50)
plt.plot(test_u, splev(test_u, spline))
plt.plot(test_u, splev(test_u, spline), 'o')
plt.cla()
plt.plot(test_u, splev(test_u, spline), 'o')
test_u = np.random.uniform(0,1,100)
plt.plot(test_u, splev(test_u, spline), 'o')
y
rand_u = np.random.uniform(0,1,1000)
plt.plot(test_u, splev(test_u, spline), 'o')
splev(test_u, spline)
plt.hist(splev(test_u, spline))
spline = splrep(bg_mean.cumsum() / bg_mean.cumsum().max(), r_unique, k=3)
splev(test_u, spline)
spline(0.5)
inv_F = interp1d(bg_mean.cumsum()/bg_mean.cumsum().max(), r_unique)
test_u
plt.cla()
plt.plot(test_u, inv_F(test_u), 'o')
bg_mean.cumsum()/bg_mean.cumsum().max()
plt.plot(test_r, splev(np.random.randint(0,1,50), spline))
splev(test_r, spline)
bg_mean.cumsum()/bg_mean.cumsum().max()
plt.plot(r_unique, bg_mean.cumsum()/bg_mean.cumsum().max(), 'o')
inv_F = interp1d(bg_mean.cumsum()/bg_mean.cumsum().max(), r_unique)
r_unique
plt.plot(bg_mean.cumsum()/bg_mean.cumsum().max(), r_unique)
bg.cumsum()
np.append(0,bg)
np.append(bg[0],bg)
np.append(bg_mean[0],bg_mean)
len(np.append(bg_mean[0],bg_mean))
len(bg_mean)
bg_mean_app = np.append(bg_mean[0],bg_mean)
plt.plot(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), r_unique)
r_unique
bg_mean
len(bg_mean)
len(bg_mean.cumsum())
plt.plot(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), r_unique)
plt.plot(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), np.append(r_unique, r_unique[-1]))
bg_mean.cumsum()/bg_mean.cumsum().max()
inv_F = interp1d(np.append(0, bg_mean.cumsum()/bg_mean.cumsum().max()), np.append(0,r_unique))
plt.hist(inv_F(test_u))
plt.cla()
plt.hist(inv_F(test_u))
test_u
test_u = np.random.uniform(0,1,1000)
plt.cla()
plt.hist(inv_F(test_u))
plt.hist(inv_F(test_u), normed=True)
plt.cla()
plt.hist(inv_F(test_u), normed=True)
plt.plot(r_unique, bg_mean/bg_mean.max())
test_u = np.random.uniform(0,1,10000)
plt.hist(inv_F(test_u), normed=True)
test_u = np.random.uniform(0,1,100000)
plt.hist(inv_F(test_u), normed=True)
plt.hist(inv_F(test_u),bins=50, normed=True)
plt.hist(inv_F(test_u),bins=48, normed=True)
plt.hist(inv_F(test_u),bins=20, normed=True)
plt.plot(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), r_unique)
plt.plot(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), np.append(r_unique, r_unique[-1]))
r_unique
bg_mean.cumsum()/bg_mean.cumsum().max()
plt.plot(r_unique, bg_mean/bg_mean.max())
plt.plot(np.append(-1, r_unique), bg_mean.cumsum()/bg_mean.cumsum().max())
plt.plot(np.append(-1, r_unique), bg_mean_app.cumsum()/bg_mean_app.cumsum().max())
plt.plot(np.append(-0.1, r_unique), bg_mean_app.cumsum()/bg_mean_app.cumsum().max())
plt.plot(np.append(0, r_unique), bg_mean_app.cumsum()/bg_mean_app.cumsum().max())
plt.plot(np.append(0, r_unique), bg_mean_app.cumsum()/bg_mean_app.cumsum().max())
inv_F = interp1d(bg_mean_app.cumsum()/bg_mean_app.cumsum().max()), np.append(0,r_unique))
inv_F = interp1d(bg_mean_app.cumsum()/bg_mean_app.cumsum().max(), np.append(0,r_unique))
plt.hist(inv_F(test_u),bins=48, normed=True)
inv_F(0)
bg_mean_app.cumsum()/bg_mean_app.cumsum().max()
bg_mean.cumsum()/bg_mean.cumsum().max()
np.append(0,bg_mean.cumsum()/bg_mean.cumsum().max())
cumsum = np.append(0,bg_mean.cumsum()/bg_mean.cumsum().max())
r_unique
plt.plot(np.append(0,r_unique), cumsum)
plt.cla()
plt.plot(np.append(0,r_unique), cumsum)
inv_F = interp1d(cumsum, np.append(0,r_unique))
inv_F(0)
inv_F(0.1)
plt.hist(inv_F(test_u),bins=48, normed=True)
test_u
plt.hist(test_u)
plt.hist(test_u, normed=True)
plt.cla()
plt.hist(test_u, normed=True)
plt.hist(inv_F(test_u),bins=48, normed=True)
plt.plot(inv_F(test_u), 'o')
plt.plot(test_u,inv_F(test_u), 'o')
plt.cla()
plt.plot(test_u,inv_F(test_u), 'o')
test_u = np.random.uniform(0,1,100)
plt.plot(test_u,inv_F(test_u), 'o')
r_unique
plt.plot(np.ones(len(r_unique))*-1, r_unique)
plt.cla()
plt.plot(test_u,inv_F(test_u), 'o')
plt.plot(np.ones(len(r_unique))*-0.25, r_unique, 'o')
plt.plot(cumsum,r_unique, 'o')
plt.plot(cumsum[1:],r_unique, 'o')
plt.plot(np.ones(len(r_unique))*-0.3, r, 'o')
plt.plot(np.ones(len(r))*-0.3, r, 'o')
test_u = np.random.uniform(cumsum.min(),1,100)
plt.plot(test_u,inv_F(test_u), 'o')
cumsum.min
cumsum.min()
cumsum = g_mean.cumsum()/bg_mean.cumsum().max()
cumsum = bg_mean.cumsum()/bg_mean.cumsum().max()
cumsum
test_u = np.random.uniform(cumsum.min(),1,100)
plt.plot(test_u,inv_F(test_u), 'o')
test_u = np.random.uniform(cumsum.min(),1,10000)
plt.plot(test_u,inv_F(test_u), 'o')
plt.hist(inv_F(test_u),bins=48, normed=True)
plt.cla()
plt.hist(inv_F(test_u),bins=48, normed=True)
plt.hist(inv_F(test_u),bins=30, normed=True)
plt.hist(inv_F(test_u),bins=15, normed=True)
plt.plot(test_u,inv_F(test_u), 'o')
plt.cla()
plt.plot(test_u,inv_F(test_u), 'o')
plt.hist(inv_F(test_u),bins=15, normed=True)
plt.hist(inv_F(test_u),bins=15)
plt.hist(inv_F(test_u),bins=14)
plt.hist(inv_F(test_u),bins=100)
plt.hist(inv_F(test_u),bins=14)
inv_F = interp1d(cumsum, r_unique, kind='cubic')
inv_F = interp1d(cumsum, r_unique, kind='cubic')
cumsum
cumsum.shape
r_unique.shape
inv_F = interp1d(cumsum, r_unique, kind='cubic')
cumsum
x = [1,2,3,4,5,6,7]
y = [7,1,3,9,3,8,7]
inv_F = interp1d(x, y, kind='cubic')
x = [1,2,3,4,5,6,7,7,7]
y = [7,1,3,9,3,8,7,8,9]
inv_F = interp1d(x, y, kind='cubic')
cumsum
r_unique
bg_mean.cumsum
bg_mean.cumsum()
bg_mean
bg_mean[-37]
bg_mean[-38]
bg_mean[:-37]
bg_mean[-38]
bg_mean_cut = bg_mean[:-37]
r_unique_cut = r_unique[:-37]
inv_F = interp1d(bg_mean_cut.cumsum()/bg_mean_cut.cumsum().max(), r_unique_cut, kind='cubic')
plt.plot(test_u, inv(test_u))
plt.plot(test_u, inv_F(test_u))
plt.cla()
plt.plot(test_u, inv_F(test_u))
plt.plot(test_u, inv_F(test_u), o'o)
plt.plot(test_u, inv_F(test_u), 'o')
bg
plt.plot(r_unique, bg_mean)
def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))
bg
bg
bg_mean
bg_men_cut
bg_mean_cut
%save try_more_sampling.py
%save 180-390 try_more_sampling.py
%history
%history -f try_more_sampling.py
def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))
from scipy.optimize import curve_fit
popt,pcov = curve_fit(gaus,r_unique_cut, bg_mean_cut, p0=[1,mean,sigma])
popt,pcov = curve_fit(gaus,r_unique_cut, bg_mean_cut, p0=[1, 0 , 5])
from scipy import asarray as ar,exp
popt,pcov = curve_fit(gaus,r_unique_cut, bg_mean_cut, p0=[1, 0 , 5])
popt
plt.plot(r_unique_cut,gaus(r_unique_cut,*popt),'ro:',label='fit')
fit_bg_gauss = gaus(r_unique_cut,*popt)
plt.plot(r_unique_cut, fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
fit_bg_gauss = gaus(np.linspace(0, r_unique_cut.max(), 100),*popt)
plt.plot(r_unique_cut, fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
plt.plot(np.linspace(0,r_unique_cut, 100), fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
plt.plot(np.linspace(0,r_unique_cut.max(), 100), fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
plt.plot(np.linspace(0,r_unique_cut, 100), fit_bg_gauss)
plt.plot(np.linspace(0,r_unique_cut.max(), 100), fit_bg_gauss)
plt.plot(np.linspace(0,r_unique_cut.max(), 100), fit_bg_gauss, 'o')
plt.plot(np.linspace(0,r_unique_cut.max(), 100), fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
fit_bg_gauss = gaus(np.linspace(0, r_unique_cut.max(), 10),*popt)
plt.plot(np.linspace(0,r_unique_cut.max(), 10), fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
fit_bg_gauss = gaus(np.linspace(0, r_unique_cut.max(), 1000),*popt)
plt.plot(np.linspace(0,r_unique_cut.max(), 1000), fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max())
inv_F = interp1d(fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max(), np.linspace(0,r_unique_cut.max(), 1000) kind='cubic')
inv_F = interp1d(fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max(), np.linspace(0,r_unique_cut.max(), 1000), kind='cubic')
test_u = np.random.uniform(fit_bg_gauss.cumsum.min(),1,10000)
test_u = np.random.uniform(fit_bg_gauss.cumsum().min(),1,10000)
plt.hist(inv_F(test_u))
plt.plot(test_u, inv_F(test_u), 'o')
test_u.min()
fit_bg_gauss.cumsum()
test_u = np.random.uniform((fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max).min(),1,10000)
cumsum = fit_bg_gauss.cumsum()/fit_bg_gauss.cumsum().max()
cumsum.min()
test_u = np.random.uniform(cumsum.min(),1,10000)
plt.plot(test_u, inv_F(test_u), 'o')
plt.hist(inv_F(test_u))
plt.hist(inv_F(test_u), normed=True)
plt.cla()
plt.hist(inv_F(test_u), normed=True)
plt.hist(inv_F(test_u),bins=30, normed=True)
plt.plot(r,bg, 'o')
plt.hist(inv_F(test_u))
plt.plot(r,bg*800, 'o')
plt.plot(r,bg, 'o')
plt.hist(inv_F(test_u), bins=30)
plt.plot(r,bg*250, 'o')
plt.cla()
plt.plot(r,bg*250, 'o')
plt.hist(inv_F(test_u), bins=30)
plt.plot(r_unique_cut,bg_mean_cut*250, 'o')
plt.cla()
plt.plot(r_unique_cut,bg_mean_cut*250, 'o')
plt.hist(inv_F(test_u), bins=30)
plt.hist(inv_F(test_u), bins=100)
plt.plot(r_unique_cut,bg_mean_cut*80, 'o')
%history -f try_more_sampling.py
