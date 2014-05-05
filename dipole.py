# -*- coding: utf-8 -*-
"""
This file gives the dipole radiation (E and B field) in the far field, the full radiation (near field + far field) and the near field radiation only

@author: manu
"""
from __future__ import division
import numpy
c=299792458.
pi=numpy.pi
mu0=4*pi*1e-7
eps0=1./(mu0*c**2)


def Hertz_dipole_ff (r, p, R, phi, f, t=0, epsr=1.):
  """
  Calculate E and B field strength of hertzian dipole(s) in the far field
  p: array of dipole moments [[px0,py0,pz0],[px1,py1,pz1],...[pxn,pyn,pzn]]
  R: array of dipole positions [[X0,Y0,Z0],[X1,Y1,Z1],...[Xn,Yn,Zn]]
  r: observation point [x,y,z]
  f: array of frequencies [f0,f1,...]
  t: time
  phi: array with dipole phase angles (0..2pi) [phi0,phi1,...,phin]
  return: fields values at observation point r at time t for every frequency in f. E and B are (3 components,number of frequencies) arrays.
  """
  nf = len(f)
  rprime = r-R  # r'=r-R
  if numpy.ndim(p) < 2:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2))
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    w = 2*pi*f  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r'x p
    rp_c_p_c_rp = numpy.cross(rprime_cross_p, rprime) # (r' x p) x r'
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[0],(nf,1))).T
    Ey = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[1],(nf,1))).T
    Ez = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[2],(nf,1))).T
    Bx = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[0],(nf,1)).T)
    By = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[1],(nf,1)).T)
    Bz = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[2],(nf,1)).T)
    E = numpy.vstack((Ex,Ey,Ez))
    B = numpy.vstack((Bx,By,Bz))
  else:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2,axis=1)) # |r'|
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    fp = numpy.tile(f,(len(magrprime),1))
    w = 2*pi*fp  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r'x p
    rp_c_p_c_rp = numpy.cross(rprime_cross_p, rprime) # (r' x p) x r'
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[:,0],(nf,1))).T
    Ey = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[:,1],(nf,1))).T
    Ez = (w**2/(c**2*magrprimep**3) * expfac)* (numpy.tile(rp_c_p_c_rp[:,2],(nf,1))).T
    Bx = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,0],(nf,1)).T)
    By = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,1],(nf,1)).T)
    Bz = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,2],(nf,1)).T)
    E = numpy.vstack((numpy.sum(Ex,axis=0),numpy.sum(Ey,axis=0),numpy.sum(Ez,axis=0)))
    B = numpy.vstack((numpy.sum(Bx,axis=0),numpy.sum(By,axis=0),numpy.sum(Bz,axis=0)))
  return E,B


def Hertz_dipole (r, p, R, phi, f, t=0, epsr=1.):
  """
  Calculate E and B field strength of hertzian dipole(s)
  p: array of dipole moments [[px0,py0,pz0],[px1,py1,pz1],...[pxn,pyn,pzn]]
  R: array of dipole positions [[X0,Y0,Z0],[X1,Y1,Z1],...[Xn,Yn,Zn]]
  r: observation point [x,y,z]
  f: array of frequencies [f0,f1,...]
  t: time
  phi: array with dipole phase angles (0..2pi) [phi0,phi1,...,phin]
  return: fields values at observation point r at time t for every frequency in f. E and B are (3 components,number of frequencies) arrays.
  """
  nf = len(f)
  rprime = r-R  # r'=r-R
  if numpy.ndim(p) < 2:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2))
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    w = 2*pi*f  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r'x p
    rp_c_p_c_rp = numpy.cross(rprime_cross_p, rprime) # (r' x p) x r'
    rprime_dot_p = numpy.sum(rprime*p)
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[0],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[0]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[0].T,(len(f),1)).T))
    Ey = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[1],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[1]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[1].T,(len(f),1)).T))
    Ez = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[2],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[2]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[2].T,(len(f),1)).T))
    Bx = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[0],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    By = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[1],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    Bz = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[2],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    E = numpy.vstack((Ex,Ey,Ez))
    B = numpy.vstack((Bx,By,Bz))
  else:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2,axis=1))
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    fp = numpy.tile(f,(len(magrprime),1))
    w = 2*pi*fp  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r' x p
    rp_c_p_c_rp = numpy.cross(rprime_cross_p, rprime) # (r' x p) x r'
    rprime_dot_p = numpy.sum(rprime*p,axis=1)
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[:,0],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,0]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,0].T,(len(f),1)).T))
    Ey = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[:,1],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,1]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,1].T,(len(f),1)).T))
    Ez = expfac*(w**2/(c**2*magrprimep**3) * (numpy.tile(rp_c_p_c_rp[:,2],(nf,1))).T+(1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,2]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,2].T,(len(f),1)).T))
    Bx = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,0],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    By = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,1],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    Bz = expfac/(magrprimep**2*c**3)*(w**2*numpy.tile(rprime_cross_p[:,2],(nf,1)).T)*(1-c/(1j*w*magrprimep))
    E = numpy.vstack((numpy.sum(Ex,axis=0),numpy.sum(Ey,axis=0),numpy.sum(Ez,axis=0)))
    B = numpy.vstack((numpy.sum(Bx,axis=0),numpy.sum(By,axis=0),numpy.sum(Bz,axis=0)))
  return E,B


def Hertz_dipole_nf (r, p, R, phi, f, t=0, epsr=1.):
  """
  Calculate E and B field strength of hertzian dipole(s)
  p: array of dipole moments [[px0,py0,pz0],[px1,py1,pz1],...[pxn,pyn,pzn]]
  R: array of dipole positions [[X0,Y0,Z0],[X1,Y1,Z1],...[Xn,Yn,Zn]]
  r: observation point [x,y,z]
  f: array of frequencies [f0,f1,...]
  t: time
  phi: array with dipole phase angles (0..2pi) [phi0,phi1,...,phin]
  return: fields values at observation point r at time t for every frequency in f. E and B are (3 components,number of frequencies) arrays.
  """
  nf = len(f)
  rprime = r-R  # r'=r-R
  if numpy.ndim(p) < 2:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2))
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    w = 2*pi*f  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r'x p
    rprime_dot_p = numpy.sum(rprime*p)
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[0]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[0].T,(len(f),1)).T))
    Ey = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[1]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[1].T,(len(f),1)).T))
    Ez = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[2]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[2].T,(len(f),1)).T))
    Bx = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[0],(nf,1)).T)*1j
    By = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[1],(nf,1)).T)*1j
    Bz = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[2],(nf,1)).T)*1j
    E = numpy.vstack((Ex,Ey,Ez))
    B = numpy.vstack((Bx,By,Bz))
  else:
    magrprime = numpy.sqrt(numpy.sum((rprime)**2,axis=1)) #|r'|
    magrprimep = numpy.tile(magrprime, (len(f),1)).T
    phip = numpy.tile(phi, (len(f),1))
    fp = numpy.tile(f,(len(magrprime),1))
    w = 2*pi*fp  # \omega
    k = w/c     # wave number
    krp = k*magrprimep  # k|r'|
    rprime_cross_p = numpy.cross(rprime, p) # r' x p
    rprime_dot_p = numpy.sum(rprime*p,axis=1) # r'.p
    expfac = numpy.exp(1j*(w*t-krp+phip.T))/(4*pi*eps0*epsr)
    Ex = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,0]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,0].T,(len(f),1)).T))
    Ey = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,1]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,1].T,(len(f),1)).T))
    Ez = expfac*((1/magrprimep**3-w*1j/(c*magrprimep**2))*(numpy.tile(3*rprime[:,2]*rprime_dot_p,(len(f),1)).T/magrprimep**2-numpy.tile(p[:,2].T,(len(f),1)).T))
    Bx = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[:,0],(nf,1)).T)*1j
    By = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[:,1],(nf,1)).T)*1j
    Bz = expfac/(magrprimep**3*c**2)*(w*numpy.tile(rprime_cross_p[:,2],(nf,1)).T)*1j
    E = numpy.vstack((numpy.sum(Ex,axis=0),numpy.sum(Ey,axis=0),numpy.sum(Ez,axis=0)))
    B = numpy.vstack((numpy.sum(Bx,axis=0),numpy.sum(By,axis=0),numpy.sum(Bz,axis=0)))
    return E,B

if __name__ == "__main__":
  from pylab import *
  #observation points
  nx=401
  xmax=2
  nz=201
  zmax=1
  x=numpy.linspace(-xmax,xmax,nx)
  y=0
  z=numpy.linspace(-zmax,zmax,nz)

  #dipole

  freq=numpy.array([1000e6])
  #dipole moment
  p=numpy.array([0,0,1e-7])
  R=numpy.array([0,0,0])
  #dipole phases
  phases_dip=0

  nt=100
  t0=1/freq/10
  t1=5/freq
  nt=int(t1/t0)
  t=numpy.linspace(t0,t1,nt)

  print("Computing the radiation...")
  fig = figure(num=1,figsize=(10,6),dpi=300)
  for k in range(nt):
    P=numpy.zeros((nx,nz))
    for i in range(nx):
      for j in range(nz):
        r=array([x[i],y,z[j]])
        E,B=Hertz_dipole (r, p, R, phases_dip, freq, t[k], epsr=1.)
        S=real(E)**2#0.5*numpy.cross(E.T,conjugate(B.T))
        P[i,j]=sqrt(sum(S))
    print('%2.1f/100'%((k+1)/nt*100))
    #Radiation diagram
    pcolor(x,z,P[:,:].T,cmap='hot')
    fname = 'img_%s' %(k)
    clim(0,5000000)
    axis('scaled')
    xlim(-xmax,xmax)
    ylim(-zmax,zmax)
    xlabel(r'$x/$m')
    ylabel(r'$z/$m')
    title(r'$t=%2.2f$ ns'%(t[k]/1e-9))
    print 'Saving frame', fname
    fig.savefig(fname+'.png',bbox='tight')
    clf()
  close()
