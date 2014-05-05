# -*- coding: utf-8 -*-
"""
Created on Thu May  1 19:15:21 2014

@author: manu
"""

from __future__ import division
from numpy import *
from numpy.random import random
from dipole import Hertz_dipole

nt=181
np=360

#theta=linspace(0,pi,nt)
theta=arccos(2*linspace(0,1,nt)-1) #unifomily separated points on theta
phi=linspace(2*pi/np,2*pi,np)
distance=100 # measurement distance in m

f0=5e6
f1=1e9
nf=200
freq=linspace(f0,f1,nf)

#random EUT
n_dip=30 #number of Herztian dipoles

#dipole positions on the EUT
a_EUT=0.5 #radius of the EUT in m
phi_posdip=2*pi*random(n_dip)
th_posdip=arccos(2*random(n_dip)-1)

R=(array([a_EUT*sin(th_posdip)*cos(phi_posdip),a_EUT*sin(th_posdip)*sin(phi_posdip),a_EUT*cos(th_posdip)])).T

#dipole moments
pmax=1e-7 #maximum dipole moment p
r_dip=pmax*random(n_dip)
phi_dip=2*pi*random(n_dip)
th_dip=arccos(2*random(n_dip)-1)
p=(array([r_dip*sin(th_dip)*cos(phi_dip),r_dip*sin(th_dip)*sin(phi_dip),r_dip*cos(th_dip)])).T
#dipole phases
phases_dip=2*pi*random(n_dip)

P=zeros((nt,np,nf))
print("Computing the radiation...")
for i in range(nt):
  for j in range(np):
    r=(array([distance*sin(theta[i])*cos(phi[j]),distance*sin(theta[i])*sin(phi[j]),distance*cos(theta[i])])).T
    E,B=Hertz_dipole (r, p, R, phases_dip, freq, t=0, epsr=1.)
    P[i,j,:]=sqrt(sum((0.5*abs(cross(E.T,(B.T))))**2,axis=1))
    #P[i,j,:]=.5*sum(abs(E)**2,axis=0) # gives the average value of the power over a period, if you want the radiated power at time t=0, please consider using 0.5*sum(real(E)**2,axis=0)
  print('%2.1f/100'%((i+1)/nt*100))

print("Printing the radiation patterns")
from pylab import *
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['legend.fontsize'] = 'medium'
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})

#Radiation diagram
fig = figure(figsize=(12, 12), dpi=50)
for f in range(nf):
  ax = fig.add_subplot(111, projection='3d', frame_on=False)
  ax._axis3don = False
  R = P[:,:,f]/(P[:,:,f]).max()
  D = P[:,:,f].max()/(P[:,:,f]).mean()
  x = R.T  * outer(cos(phi), sin(theta))
  y = R.T  * outer(sin(phi), sin(theta))
  z = R.T  * outer(ones_like(phi), cos(theta))
  ax.plot_surface(x, y, z,  rstride=2, cstride=2,color='w',\
      linewidth=0.6,shade=False)
  max_radius = 0.7
  print('f= %2.1f MHz' %(freq[f]/1e6))
  title(r'$f= %2.1f$ MHz, $D \approx %2.1f$' %(freq[f]/1e6,D),fontsize=20)
  for axis in 'xyz':
      getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
  fname = 'rp_%s' %(f)
  print 'Saving frame', fname
  fig.savefig(fname+'.png',bbox='tight')
  #fig.savefig(fname+'.svg',bbox='tight')
  #fig.savefig(fname+'.pdf',bbox='tight')
  clf()
close()
