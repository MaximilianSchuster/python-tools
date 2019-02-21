# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 10:33:05 2019

@author: schuster
"""
#%%
import numpy as np
import os
#import configparser
import scipy
from scipy import io
from scipy import signal
import insect_tools as insect
import wabbit_tools as w2p

# import the geometric variables from the ini-file

oF_path1 = '/home/schuster/Masterarbeit/Matlab/OneFlame/Q_short.mat' # sfb-pdc geometry

oF_path2 = '/home/schuster/Masterarbeit/Matlab/OneFlame/Q_short2.mat' # const geometry
daten2 = scipy.io.loadmat(oF_path1)
q = daten2['Q_short']

rho = q[:,:,0]
u = q[:,:,1]/rho
p = q[:,:,2]

rho_OF= rho[:,1039:3999]
u_OF=u[:,1039:3999]
p_OF=p[:,1039:3999]

rho_OF_comp = np.transpose(rho_OF[:,0:2960:40])
u_OF_comp = np.transpose(u_OF[:,0:2960:40])
p_OF_comp = np.transpose(p_OF[:,0:2960:40])
wabbit_path1 = '/home/schuster/WABBIT/flow_var_pdc'
wabbit_path1 = '/work/schuster/wabbit_plenum_new/flow_var/'
wabbit_path2 = '/home/schuster/WABBIT/flow_var_const'

rho_WABBIT = w2p.space_time_slice(wabbit_path1,'rho', 259, None)[0:74,:]
u_WABBIT = w2p.space_time_slice(wabbit_path1,'Ux', 259, None)[0:74,:]
p_WABBIT = w2p.space_time_slice(wabbit_path1,'p', 259, None)[0:74,:]

#%%
x_lim_4 = 512
x_lim_3 = 256
c_lo_rho = 0
c_hi_rho = 7
c_lo_p = 11
c_hi_p = 16
c_lo_u = 0
c_hi_u = 1300
import matplotlib.pyplot as plt
plt.figure()
plt.axis('off')

plt.subplot(321)
plt.pcolor(rho_OF_comp)
plt.xlim([0,x_lim_3])
plt.ylim([0, 69])
plt.title('rho OneFlame')
plt.colorbar()
plt.clim([c_lo_rho, c_hi_rho])

plt.subplot(323)
plt.pcolor(np.log(p_OF_comp))
plt.xlim([0,x_lim_3])
plt.ylim([0, 69])
plt.title('P OneFlame')
plt.colorbar()
plt.clim([c_lo_p, c_hi_p])


plt.subplot(325)
plt.pcolor(u_OF_comp)
plt.xlim([0,x_lim_3])
plt.ylim([0, 69])
plt.title('Ux OneFlame')
plt.colorbar()
plt.clim([c_lo_u, c_hi_u])


plt.subplot(322)
plt.pcolor(rho_WABBIT)
plt.xlim([0,x_lim_4])
plt.ylim([0, 69])
plt.title('rho WABBIT')
plt.colorbar()
plt.clim([c_lo_rho, c_hi_rho])


plt.subplot(324)
plt.pcolor(np.log(p_WABBIT))
plt.xlim([0,x_lim_4])
plt.ylim([0, 69])
plt.title('P WABBIT')
plt.colorbar()
plt.clim([c_lo_p, c_hi_p])


plt.subplot(326)
plt.pcolor(u_WABBIT)
plt.xlim([0,x_lim_4])
plt.ylim([0, 69])
plt.colorbar()
plt.title('Ux WABBIT')
plt.clim([c_lo_u, c_hi_u])


plt.tight_layout()
plt.show()
plt.savefig('WABBIT_OneFlame.pdf')




#
##%% ARE THE INI CONDS RIGHT COMP
#
#import wabbit_tools as w2p
#
#pr = w2p.read_wabbit_hdf5('/home/schuster/WABBIT/p_000000000000.h5')
#rho = w2p.read_wabbit_hdf5('/home/schuster/WABBIT/rho_000000000000.h5')
#vel = w2p.read_wabbit_hdf5('/home/schuster/WABBIT/Ux_000000000000.h5')
#k = w2p.dense_matrix(rho[1], rho[2], rho[4], rho[5], dim=2)
#k = k[0]
##%%
#plt.figure()
#plt.plot(k[128,:])
#plt.plot(rho_OF_comp[0,:])
#plt.show()



#%%
#rho = w2p.read_wabbit_hdf5('/home/schuster/WABBIT/flow_var/rho_000000000150.h5')
#k = w2p.dense_matrix(rho[1], rho[2], rho[4], rho[5], dim=2)
#k = k[0]


##simple line plot comp.
