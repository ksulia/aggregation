import matplotlib.pyplot as plt      #import 2D plotting library
import matplotlib as m
from matplotlib import colors, ticker, cm, mlab
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import numpy as np                   #import scientific computing library
import pylab as pl
import os                            #import operating system library
from netCDF4 import Dataset          #import dataset library from netCDF pkg
import glob
import sys
from pylab import *
import matplotlib.font_manager as font_manager
import warnings
warnings.filterwarnings("ignore")
np.set_printoptions(threshold=np.inf)


infile = 'OUTPUT.NC'
fh = Dataset(infile,mode = 'r')

time = (fh.variables['TIMES'][:])
qi = (fh.variables['QICE'][:])
qc = (fh.variables['QCLOUD'][:])
qv = (fh.variables['QVAPOR'][:])
qs = (fh.variables['QSNOW'][:])
qr = (fh.variables['QRAIN'][:])
ni = (fh.variables['QNICE'][:])
ns = (fh.variables['QNSNOW'][:])
nr = (fh.variables['QNRAIN'][:])
nc = (fh.variables['QNCLOUD'][:])
aice = (fh.variables['AICE'][:])
cice = (fh.variables['CICE'][:])
th = (fh.variables['THETA'][:])
p = (fh.variables['PRESS'][:])
relh = (fh.variables['RELH'][:])
phi = (fh.variables['PHI'][:])
rhoice =  (fh.variables['RHOICE'][:])
nuc = (fh.variables['ICENUC'][:])
dep = (fh.variables['ICEDEP'][:])
sub = (fh.variables['ICESUB'][:])

temp = th*(p/1.e5)**(287.05/1005)
numt = time.shape[0]
ai = np.zeros((numt))
ci = np.zeros((numt))
wh = np.where((aice > 0.0) & (cice > 0.) & (ni > 0.0))
ai[wh] = (aice[wh]*aice[wh]/(cice[wh]*ni[wh]*4.))**(1./3.)
ci[wh] = (cice[wh]*cice[wh]/(aice[wh]*ni[wh]*4.))**(1./3.)


gs = gridspec.GridSpec(5,3,width_ratios=[1,1,1],hspace=0.7)
axis_font = {'size':'14'}
TITLE = str(int(temp[1]-273.15))+' $^o$C'


ax1 = plt.subplot(gs[0])
P1 = ax1.plot(time,relh)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('RH', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')

ax1 = plt.subplot(gs[1])
P1 = ax1.plot(time,temp-273.15)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('Temp (C)', **axis_font)

ax1 = plt.subplot(gs[2])
P1 = ax1.plot(time,phi)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('Aspect Ratio', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
plt.title(TITLE)

ax1 = plt.subplot(gs[3])
P1 = ax1.plot(time,qv*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('q$_{vapor}$ g kg$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qv) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[4])
P1 = ax1.plot(time,qi*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('q$_{ice}$ g kg$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qi) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[5])
P1 = ax1.plot(time,qc*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('q$_{cloud}$ g kg$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qc) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[6])
P1 = ax1.plot(time,qs*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('q$_{snow}$ g kg$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qs) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[7])
P1 = ax1.plot(time,qr*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('q$_{rain}$ g kg$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qr) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[8])
P1 = ax1.plot(time,ai*1.e6)
P2 = ax1.plot(time,ci*1.e6)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel(r'axis length $\mu$m', **axis_font)
plt.legend(['a','c'])
ax1.ticklabel_format(style='sci', scilimits=(-4,4), axis='y')
if any(ai) > 0: ax1.set_yscale('log')


#ax1 = plt.subplot(gs[9])
#P1 = ax1.plot(time,ci*1.e6)
#plt.xlabel('Time (s)', **axis_font)
#plt.ylabel('c$_{ice}$ g kg$^{-1}$', **axis_font)
#ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
#if any(ci) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[9])
P1 = ax1.plot(time,dep*1000.)
P2 = ax1.plot(time,abs(sub)*1000.)
P3 = ax1.plot(time,nuc*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('Ice Souce', **axis_font)
plt.legend(['Dep','Sub','Nuc'])
ax1.ticklabel_format(style='sci', scilimits=(-4,4), axis='y')
if any(dep) > 0 or any(abs(sub)) > 0 or any(nuc) > 0: ax1.set_yscale('log')

ax1 = plt.subplot(gs[10])
P1 = ax1.plot(time,rhoice)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel(r'$\rho_{ice}$  kg m$^{-3}$', **axis_font)

ax1 = plt.subplot(gs[11])
P1 = ax1.plot(time,ni*1000.)
plt.xlabel('Time (s)', **axis_font)
plt.ylabel('n$_{ice}$ L$^{-1}$', **axis_font)
ax1.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
if any(qi) > 0: ax1.set_yscale('log')


plt.show()
fh.close()
