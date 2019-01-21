#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 16:44:23 2019

@author: max
"""

### THIS CODE CONTAINS LOCAL PATHS, BEFORE USING IT ON ANOTHER MACHINE CHANGE THE PATHS or type in geometric variables manually

import numpy as np
import os
import configparser
import scipy
from scipy import io
from scipy import signal
import insect_tools as insect
# import the geometric variables from the ini-file


config = configparser.ConfigParser()
config.read("/home/max/Masterarbeit/Python/Code Dev/ion_funnel_light.ini")

#-------------------------READ PARAMETERS FROM INI FILE------------------------------------
#------------------------------------------------------------------------------------------ 
# case selection
# 1 = first try with post flame distribution
# 2 = DDT setup
# 3 = acoustic pressure pulse
case = 2;
#GEOMETRY
Lx=config.get("Domain", "domain_size")
Lx=float(Lx[0:3])
Lx = 0.4
Ly=config.get("Domain", "domain_size")
Ly=float(Ly[4:7])
Ly = 0.2
D_pip = config.get("plenum", "diameter_pip")
D_pip = D_pip[:-1]
D_pip =float(D_pip)
D_ple= config.get("plenum", "diameter_ple")
D_ple= float(D_ple[:-1])
D_ple = 0.2
R_pip=D_pip/2
R_ple=D_ple/2
L_pip = float(config.get("plenum","length_pip")[:-1])
L_pip = 0.4
L_ple = float(config.get("plenum","length_ple")[:-1])
wall_thickness= config.get("plenum", "wall_thickness")
wall_thickness = float(wall_thickness[:-1])*Lx
wall_thickness = 0.00625
### INITIAL VALUES
rho0=1.204
p0= 100000
C_eta = 1e-7
C_sp = 1e-7
# READ FIELD VALUES FROM MATLAB MATRIX
daten=scipy.io.loadmat('/home/max/Masterarbeit/Python/Code Dev/data.mat')
daten2 = scipy.io.loadmat('/home/max/Masterarbeit/Matlab/Rho_u_p_Y.mat')
daten3 = np.asarray(daten2['Rho_u_p_Y'])
p_t1 = np.zeros((192,1))
u_t1 = np.zeros((192,1))
rho_t1 = np.zeros((192,1))
zz = 0;
for i in range(0,3072,16):
    rho_t1[zz,0]=daten3[i,0]
    u_t1[zz,0]=daten3[i,1]
    p_t1[zz,0]=daten3[i,2]
    zz = zz+1




daten = daten['data']
p_t0 = daten[:,0]
u_t0 = daten[:,1]
rho_t0 = daten[:,2]
# WABBIT MESH PARAMETERS
Bs = int(config.get("Blocks", "number_block_nodes")[:2])
#Bs = 17 # blocksize
lvl = int(config.get("Blocks", "min_treelevel")[:-1])
g = int(config.get("Blocks", "number_ghost_nodes")[:-1])
#------------------------FILL P/RHO/T/mask VALUES IN MATRICES------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------


Ds = Bs*2**lvl - 2**lvl+1
mask = np.asarray(np.zeros((Ds,Ds)))
P = np.asarray(p0*np.ones((Ds,Ds)))
rho = np.asarray(rho0*np.ones((Ds,Ds)))
u =np.asarray(np.zeros((Ds,Ds)))
v = np.asarray(np.zeros((Ds,Ds)))
pip_width=0    
pip_length=0
# SET UP SPONGE AND COUNT DIMENSION OF INICOND MATRIX
for iy in range(0, Ds):
    y = Ly/Ds*iy
    r = abs(Ly/2-y)
    #print(r)
    #print(y)
    #y = (iy-(g+1)) * dx(2) + x0(2)
    #r = abs(y-domain_size(2)*0.5)
    if r < R_pip:
        pip_width=pip_width+1;
    
    for ix in range(0, Ds):
            #x = (ix-(g+1)) * dx(1) + x0(1)
            x = Lx/Ds*ix
            #print(x)
            #r = abs(Ly/2-y)
            #print(r)
            # Plenum Bounds
            if x<=L_pip and r>= R_pip:
                mask[iy,ix]=1
            if x<=wall_thickness:
                mask[iy,ix]=1
            # quadratic power law determines strength of the sponge in radial plenum direction
            if x>L_pip and abs(r)>Ly/2.0-wall_thickness:
                x_min=Ly/2.0-wall_thickness
                x_max=Ly/2.0
                mask[iy,ix]=((abs(r)-x_min)/(x_max-x_min))**2
            # quadratic power law determines strength of the sponge in RIGHT BOUNDARY DOMAIN (OUTLET)
            if x>Lx-wall_thickness and abs(r)/2.0-wall_thickness: 
                x_min=Lx-wall_thickness
                x_max=Lx
                mask[iy,ix]=((x-x_min)/(x_max-x_min))**2
            # multiplying power law for "bottom right corner" and "top right corner" of the domain (--> x^4 Law?)
            if x>Lx-wall_thickness and abs(r)>=Ly/2.0-wall_thickness:
                x_min=Ly/2.0-wall_thickness
                x_max=Ly/2.0
                P1=((abs(r)-x_min)/(x_max-x_min))**2
                x_min=Lx-wall_thickness
                x_max=Lx
                P2=((x-x_min)/(x_max-x_min))**2
                mask[iy,ix] = P1*P2
pip_start=0
for ix in range(0, Ds):
     x = Lx/Ds*ix

     
     if x < L_pip and x>wall_thickness:
         pip_length=pip_length+1
     if x < wall_thickness:
         pip_start = pip_start+1
#            if x < L_pip and iy == (Ds+1)/2    :
#                P[iy,ix] = 10
for ix in range(0, Ds):
    for iy in range(0,Ds):
        if mask[ix,iy]==1:
          rho[ix,iy] = rho0
          P[ix,iy] = p0
          u[ix,iy] = 0
          v[ix,iy] = 0
          #mask[ix,iy]= mask[ix,iy]/C_sp
        elif mask [ix,iy] > 0 and mask[ix,iy] < 1:
            rho[ix,iy] = rho0
            P[ix,iy] = p0
            u[ix,iy] = 0
            v[ix,iy] = 0
           # mask[ix,iy] = mask[ix,iy]/C_sp

#pip_length=pip_length-pip_start          
ad_flame_p = np.asarray(np.zeros((pip_width,pip_length)))
ad_flame_rho = np.asarray(np.zeros((pip_width,pip_length)))
ad_flame_u = np.asarray(np.zeros((pip_width,pip_length)))
ad_flame_v = np.asarray(np.zeros((pip_width,pip_length)))
#ad_flame[int(pip_width/2),:] = 100*np.random.rand(1,pip_length)
#ad_flame[int(pip_width/2)+1,:]=ad_flame[int(pip_width/2),:]
yy = np.linspace(-1, 1, int(pip_width))
cc = np.cosh(yy)
mid = int(pip_width/2)
rr = np.around(np.linspace(0,191,pip_length))
for iy in range(int(pip_length)):#(0,int(pip_length-1)):
    kk = int(rr[iy])
    ### Fill in fully evolved profiles
    if case == 1: 
        ad_flame_p[:,iy] =np.ones(int(pip_width))*p_t0[kk]#*p_t0[k]#+ np.cosh(k[0]*np.ones(1,pip_length))
        ad_flame_rho[:,iy]=np.ones(int(pip_width))*rho_t0[kk]
        ad_flame_u[:,iy]=(-cc+cc[0])/(-cc[mid]+cc[0])*u_t0[kk]
    ### Fill in semi interrupted Y-Distribution
    elif case==2:
        ad_flame_p[:,iy] =np.ones(int(pip_width))*p_t1[kk]#*p_t0[k]#+ np.cosh(k[0]*np.ones(1,pip_length))
        ad_flame_rho[:,iy]=np.ones(int(pip_width))*rho_t1[kk]
        ad_flame_u[:,iy]=(-cc+cc[0])/(-cc[mid]+cc[0])*u_t1[kk]

        
        
        
nx_grid_min=pip_start
nx_grid_max=pip_length+pip_start
ny_grid_min=(Ds+1)/2-1-pip_width/2
ny_grid_max=(Ds+1)/2-1+pip_width/2

P[int(ny_grid_min):int(ny_grid_max),int(nx_grid_min):int(nx_grid_max)] = ad_flame_p  
rho[int(ny_grid_min):int(ny_grid_max),int(nx_grid_min):int(nx_grid_max)] = ad_flame_rho  
u[int(ny_grid_min):int(ny_grid_max),int(nx_grid_min):int(nx_grid_max)] = ad_flame_u
#%%
if case ==3:
        rho0 =1.225;
        p00 = 101325;
        puls = signal.gaussian(19,10)*150;
        Ppuls = np.ones([19,19]);
        P = np.ones([257,257])*p00;
        puls = p00+puls;
        P[118:137,118:137]=puls;
        rhopuls = np.power(puls/p00,1/1.4)*rho0;
        rho = np.ones([257,257])*rho0;
        rho[118:137,118:137]=rhopuls;
        u = np.zeros([257,257]);
        v = u;



NdF = 4 # number of datafields

dx = 2**(-lvl) * Lx / (Bs-1) # spacing
dy = 2**(-lvl)*Ly/(Bs-1)

U = np.transpose(u)
V = np.transpose(v)
RHO = np.transpose(rho)
P = np.transpose(P)

Q1 = np.zeros([1,256,256])
Q1[:,0:256,0:256]=P[0:256,0:256]
Q11 = insect.write_flusi_HDF5('p_00.h5', 0, [0, 0.4,0.2], Q1, viscosity=0.0, origin=np.array([0.0,0.0,0.0]))

Q2 = np.zeros([1,256,256])
Q2[:,0:256,0:256]=U[0:256,0:256]
Q2 = insect.write_flusi_HDF5('Ux_00.h5', 0, [0, 0.4,0.2], Q2, viscosity=0.0, origin=np.array([0.0,0.0,0.0]))

Q3 = np.zeros([1,256,256])
Q3[:,0:256,0:256]=V[0:256,0:256]
Q33 = insect.write_flusi_HDF5('Uy_00.h5', 0, [0, 0.4,0.2], Q3, viscosity=0.0, origin=np.array([0.0,0.0,0.0]))

Q4 = np.zeros([1,256,256])
Q4[:,0:256,0:256]=RHO[0:256,0:256]
Q44 = insect.write_flusi_HDF5('rho_00.h5', 0, [0, 0.4,0.2], Q4, viscosity=0.0, origin=np.array([0.0,0.0,0.0]))

#%% CONVERT FLUSI TO WABBIT 
os.system("/home/max/FORK_SCHUSTER/WABBIT/wabbit-post --flusi-to-wabbit p_00.h5 p_0.h5 33")
os.system("/home/max/FORK_SCHUSTER/WABBIT/wabbit-post --flusi-to-wabbit rho_00.h5 rho_0.h5 33")
os.system("/home/max/FORK_SCHUSTER/WABBIT/wabbit-post --flusi-to-wabbit Ux_00.h5 Ux_0.h5 33")
os.system("/home/max/FORK_SCHUSTER/WABBIT/wabbit-post --flusi-to-wabbit Uy_00.h5 Uy_0.h5 33")