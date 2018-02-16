import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('ggplot')

import astropy.table
from astropy.table import Table
tugreen = '#73ad14'




def animate_trigger_citerion_Comparison(path_to_raw_file,path_to_denoised_file,Output):

    cube_table = Table.read(path_to_raw_file, path='data')
    denoised_table = Table.read(path_to_denoised_file, path='data')
    Cube = cube_table['cube'][0]
    Den = denoised_table['cube_smoothed'][0]

    Trigger_Crit = np.zeros(60)
    Trigger_Crit_Den = np.zeros(60)
    x = np.arange(0,60)
    for i in x:
        Trigger_Crit[i] = Cube[i].max().max()
        Trigger_Crit_Den[i] = Den[i].max().max()

    y_Cube = Trigger_Crit
    y_Den = Trigger_Crit_Den

    fig, ax = plt.subplots(2,sharex=True)
    line1, = ax[0].plot(x, y_Cube, color='crimson', label='Raw_Cube')
    ax[0].set_ylabel('trigger_criterion')
    ax[0].legend()

    line2, = ax[1].plot(x, y_Den, color=tugreen, label='Denoised_Cube')
    ax[1].set_xlabel('slices')
    ax[1].set_ylabel('trigger_criterion')
    ax[1].legend()

    def update(num, x, y_Den,y_Cube, line1, line2):
        line1.set_data(x[:num], y_Cube[:num])
        line1.axes.axis([0, 60, -1, y_Cube.max()+10])
        line2.set_data(x[:num], y_Den[:num])
        line2.axes.axis([0, 60, -1, y_Cube.max()+10])
        return line1,line2


    ani = animation.FuncAnimation(fig, update, np.arange(0,60), fargs=[x, y_Den,y_Cube, line1,line2],
                              interval=100, blit = False, repeat=False)

    ani.save('%s.gif'%(Output))

animate_trigger_citerion_Comparison('../build/n4_s60_p2_cube.hdf5',
                                        '../build/n4_s60_p2_denoised.hdf5',
                                        '../Wavelet_Gauss')
