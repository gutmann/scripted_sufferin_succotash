#!/usr/bin/python
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.

#
# 2D flow around a cylinder
#

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt; from matplotlib import cm

from numba import jit
import time

import xarray as xr
ds=xr.open_dataset("profile.nc")
topography = ds.data.values[300:750]
topography -= topography.min()

nx = int(topography.shape[0])
ny = int(topography.max()+200)

###### Flow definition #########################################################
maxIter = 200000 # Total number of time iterations.
Re      = 200.0  # Reynolds number.
# nx = 520; ny = 180;
ly=ny-1.0; q = 9 # Lattice dimensions and populations.
uLB     = 0.04                       # Velocity in lattice units.
r=20
nulb    = uLB*r/Re; omega = 1.0 / (3.*nulb+0.5); # Relaxation parameter.
print(omega)
###### Lattice Constants #######################################################
c = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
t = 1./36. * np.ones(q)                                   # Lattice weights.
t[np.asarray([norm(ci)<1.1 for ci in c])] = 1./9.; t[0] = 4./9.
noslip = np.array([c.tolist().index((-c[i]).tolist()) for i in range(q)])
i1 = np.arange(q)[np.asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
i2 = np.arange(q)[np.asarray([ci[0]==0 for ci in c])] # Vertical middle.
i3 = np.arange(q)[np.asarray([ci[0]>0  for ci in c])] # Unknown on left wall.
i4 = np.arange(q)[np.asarray([ci[1]<0  for ci in c])] # Unknown on bottom wall.
i5 = np.arange(q)[np.asarray([ci[1]==0 for ci in c])] # Horizontal middle
i6 = np.arange(q)[np.asarray([ci[1]>0  for ci in c])] # Unknown on top wall.

###### Function Definitions ####################################################
# sumpop = lambda fin: sum(fin,axis=0) # Helper function for density computation.

###### Setup: obstacle and velocity inlet with perturbation ########
# obstacle = np.empty((nx,ny), dtype=bool)
# obstacle[:] = False
# obstacle[100,100:] = True
# obstacle[200,100:] = True
# obstacle[100:200,100] = True
#
# obstacle[100:200,100:] = True
# obstacle[:,-1] = True # set a no slip bottom boundary

@jit
def define_obstacle(obstacle, topography):
    for i in range(obstacle.shape[0]):
        for j in range(obstacle.shape[1]):
            obstacle[i,j] = j<=topography[i]

obstacle = np.empty((nx,ny), dtype=bool)
define_obstacle(obstacle, topography)

# obstacle = np.array(np.fromfunction(lambda x,y: y <= topography[np.int(x)], (nx,ny)))
obstacle[:,-1] = True # set a no slip bottom boundary


# Create a cylinder obstacle
# Coordinates of the cylinder.
# cx = nx/4; cy=ny/2; r=ny/9;
# obstacle = np.array(np.fromfunction(lambda x,y: (x-cx)**2+(y-cy)**2<r**2, (nx,ny)))

vel = np.fromfunction(lambda d,x,y: (1-d)*uLB*(1.0+1e-4*np.sin(y/ly*2*np.pi)),(2,nx,ny))


xpoints=[]
ypoints=[]
for x in range(nx):
    for y in range(ny):
        if obstacle[x,y]:
            xpoints.append(x)
            ypoints.append(y)

xpoints = np.array(xpoints).astype(np.int)
ypoints = np.array(ypoints).astype(np.int)


@jit(nopython=True)
def equilibrium(rho,u, feq, cu):              # Equilibrium distribution function.

    # cu   = 3.0 * np.dot(c,np.transpose(u,(1,0,2)))
    # cu = np.zeros((c.shape[0],u.shape[1],u.shape[2]))
    for k in range(c.shape[0]):
        for i in range(u.shape[1]):
            for j in range(u.shape[2]):

                # if not obstacle[i,j]:
                cu[k,i,j] = c[k,0] * u[0,i,j]
                cu[k,i,j] += c[k,1] * u[1,i,j]
    cu *= 3

    usqr = 3./2.*(u[0]**2+u[1]**2)

    for i in range(q):
        for x in range(u.shape[1]):
            for y in range(u.shape[2]):
                # if not obstacle[x,y]:
                feq[i,x,y] = rho[x,y]*t[i]*(1.+cu[i,x,y]+0.5*cu[i,x,y]**2-usqr[x,y])

    # return feq


@jit#(nopython=True)
def update(fin, u, cu, feq, xpts, ypts, noslip, obstacle):
    fin[i1,-1,:] = fin[i1,-2,:]     # Right wall: outflow condition.
    rho = np.sum(fin,axis=0)
    # rho = fin[0]
    # for i in range(1,q):
    #     rho += fin[i]            # Calculate macroscopic density and velocity.

    # u = np.dot(c.transpose(), fin.transpose((1,0,2)))/rho
    for i in range(2):
        for j in range(nx):
            for k in range(ny):
                u[i,j,k] = 0
                # if not obstacle[j,k]:
                for n in range(q):
                    u[i,j,k] += (c[n,i] * fin[n,j,k])
                u[i,j,k] /= rho[j,k]

    u[:,0,:] =vel[:,0,:] # Left wall: compute density from known populations.
    rho[0,:] = 0
    for i in range(len(i2)):
        rho[0,:] += 1./(1.-u[0,0,:]) * (fin[i2[i],0,:] + 2. * fin[i1[i],0,:])

    equilibrium(rho,u, feq, cu)#, obstacle) # Left wall: Zou/He boundary condition.
    fin[i3,0,:] = fin[i1,0,:] + feq[i3,0,:] - fin[i1,0,:]
    fout = fin - omega * (fin - feq)  # Collision step.


    for i in range(q):
        # fout[i,obstacle] = fin[noslip[i],obstacle]
        for j in range(len(xpoints)):
            fout[i,xpts[j],ypts[j]] = fin[noslip[i],xpts[j],ypts[j]]

    # for i in range(q): # Streaming step.
    #     fin[i,:,:] = np.roll(
    #                         np.roll(
    #                             fout[i,:,:],c[i,0],axis=0), c[i,1],axis=1)
    for i in range(q): # Streaming step.
        for x in range(nx):
            for y in range(ny):
                # if not obstacle[x,y]:
                fin[i,x,y] = fout[i,np.mod(x-c[i,0],nx),np.mod(y-c[i,1],ny)]

    # return u


def main():

    feq = np.zeros((q,nx,ny))
    cu = np.zeros((q,nx,ny))
    u = np.zeros((2,nx,ny))
    equilibrium(np.ones((nx,ny)),np.array(vel), feq, cu)#, obstacle)
    fin = feq.copy()

    xcoord = np.arange(x)
    ycoord = np.arange(y)
    xcoord, ycoord = np.meshgrid(xcoord, ycoord)

    ###### Main time loop ##########################################################
    plt.clf()

    # plt.subplot(4,1,1)
    # img1 = plt.imshow(feq.sum(axis=0).transpose()-0.999,cmap=cm.RdBu)
    # plt.colorbar()
    # plt.clim(-0.1,0.1)
    #
    # plt.subplot(4,1,2)
    # img2 = plt.imshow(feq.sum(axis=0).transpose()-0.999,cmap=cm.RdBu)
    # plt.colorbar()
    # plt.clim(-0.1,0.1)
    #
    # plt.subplot(4,1,3)
    # img3 = plt.imshow(feq.sum(axis=0).transpose()-0.999,cmap=cm.Reds)
    # plt.colorbar()
    # plt.clim(0,0.13)
    #
    # ax=plt.subplot(4,1,4)
    ax=plt.subplot(1,1,1)


    s=2
    steps_between_plots = 50
    t1=time.time()
    for t in range(maxIter):
        update(fin, u, cu, feq, xpoints, ypoints, noslip, obstacle)

        if (t%steps_between_plots==0): # Visualization
            t2=time.time()
            print(t2-t1)
            # plt.clf()
            # img1.set_data(u[0].transpose())
            # img2.set_data(u[1].transpose())
            # img3.set_data(np.sqrt(u[0]**2+u[1]**2).transpose())
            ax.clear()
            ax.imshow(np.sqrt(u[0]**2+u[1]**2).transpose(),cmap=cm.Reds, origin='lower')
            # ax.colorbar()
            # plt.clim(0,0.3)
            u[0][obstacle]=0
            u[1][obstacle]=0
            # u[:,100:200,100:]=0
            if s>1:
                ax.quiver(xcoord[::s,::s],ycoord[::s,::-s],u[0,::s,::s].transpose(),u[1,::s,::s].transpose(),width=0.0005,scale=8)
            else:
                ax.quiver(u[0].transpose(),u[1].transpose(),width=0.0003,scale=10)

            plt.ylim(50,250)
            # plt.xlim(300,650)
            plt.plot(topography, color="black")
            plt.pause(0.01)
            plt.savefig("output/vel3.{:04}.png".format(int(t/steps_between_plots)))
            t1=time.time()

if __name__ == '__main__':
    main()
