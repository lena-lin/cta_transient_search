from astropy.io import fits
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt

irf = fits.open('/home/lena/Dokumente/CTA/caldb/data/cta/prod3b/bcf/South_z20_average_100s/irf_file.fits')
bg = irf['BACKGROUND'].data['BGD'][0].sum(axis=0)
bg = bg.flatten()

x = (irf['BACKGROUND'].data['DETX_HI'] + irf['BACKGROUND'].data['DETX_LO'])/2
y = (irf['BACKGROUND'].data['DETY_HI'] + irf['BACKGROUND'].data['DETY_LO'])/2

X, Y = np.meshgrid(x, y)
r = np.sqrt(X ** 2 + Y ** 2)
r = r.flatten()
sort_idx = np.argsort(r)
r = r[sort_idx]
bg = bg[sort_idx]

bg_cut = bg[:np.where(bg == 0)[0].min()]
r_cut = r[:np.where(bg == 0)[0].min()]

df_bg = pd.DataFrame(data={'r': r_cut, 'bg': bg_cut})
df_bg_grouped = df_bg.groupby(['r']).mean()

r_unique = df_bg_grouped.index.values
bg_mean = df_bg_grouped['bg'].values

bg_mean = np.append(bg_mean[0], bg_mean)
r_unique = np.append(0, r_unique)

r_even = np.linspace(r_unique.min(), r_unique.max(), 100)
spline_f = splrep(r_unique, r_unique * bg_mean, xb=0, xe=6, k=3)
cumsum = splev(r_even, spline_f).cumsum()/splev(r_even, spline_f).cumsum().max()
inv_F = interp1d(cumsum, r_even)

u = np.random.uniform(cumsum.min(), 1, 500)
r_sampled = inv_F(u)

phi = np.random.uniform(0, 2 * np.pi, len(u))

ra = r_sampled * np.cos(phi)
dec = r_sampled * np.sin(phi)

plt.figure()
plt.hist2d(r_sampled * np.cos(phi), r_sampled * np.sin(phi), bins=[36, 36])

hist = np.histogram2d(ra, dec, bins=[36, 36])
bg_sampled = hist[0].flatten()
x = (hist[1][:-1] + hist[1][1:])/2
y = (hist[2][:-1] + hist[2][1:])/2
X, Y = np.meshgrid(x, y)
r_sampled = np.sqrt(X ** 2 + Y ** 2)
r_sampled = r_sampled.flatten()
sortx = np.argsort(r_sampled)
r_sampled = r_sampled[sortx]
bg_sampled = bg_sampled[sortx]
