import numpy as np
import matplotlib.pyplot as plt
import os
import f90nml
import analytical_solution as aso
import signal_processing as sp
import directory_edit as de
import phys_cons as pc

de.make_dir("fig")
params = f90nml.read("param.inp")
nx   =  params["space"]["nxx"]
dx   =  params["space"]["dx"]
io   =  params["output"]["io"]
dy   =  dx
c    =  pc.c
z0   =  pc.z0
dt   =  params["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
t0   =  params["scatt"]["tau0"]
prad    =  params["plasma"]["prad"]
wp = params["plasma"]["wp"]
fp = wp/2/np.pi
nu = params["plasma"]["nu"]
d = 2*prad*dx
f1   =  c/(10*dx)
alpha=  16/t0/t0
amp  =  params["scatt"]["amp"]
N = params["time"]["nstep"]
N = N + 1

t = np.arange(0,N*dt,dt)
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
k0  = 2*np.pi*freq/c

iez = sp.input_signal("eyin.txt")
iez = sp.add_zero(iez,N)
fiez = sp.fft(iez)
sp.plot_signal(t,iez,"incident wave",xmax=0.5e-9)
sp.show(draw=False,save=False,save_name="incident_wave")
sp.plot_fsignal(freq,fiez,xmin=1e6,xmax=80e9)
sp.show(draw=False,save=False,save_name="fincident_wave")

rez = sp.input_signal("eyinc.txt")
frez = sp.fft(rez)
sp.plot_signal(t,rez,xmax=0.5e-9,title="scatterd wave")
sp.show(draw=True,save=True,save_name="scattered_wave")
sp.plot_fsignal(freq,frez,xmin=1e6,xmax=80e9,title="fscatterd_wave")
sp.show(draw=True,save=False,save_name="fscatterd_wave")

rcs = sp.rcs(fiez,frez)
afreq,ref,trs = aso.slab_mag_perpendicular(wp,3e11,nu,d)
sp.plot_fsignal(freq,rcs,color="red",marker="o",xmin=1e9,xmax=80e9,ymin=-45,ymax=0,label="simulation")
if(io<int(nx/2)):
    sp.plot_fsignal(afreq,ref,color="blue",ylabel="Reflection Coefficient Magnitude(dB)",label="analytical")
    plt.vlines(fp,ymin=-45,ymax=0,linestyle="dashed",color="black")
    sp.show(save=True,save_name="reflection")
else:
    sp.plot_fsignal(afreq,trs,color="blue",ylabel="Transmission Coefficient Magnitude(dB)",label="analytical")
    plt.vlines(fp,ymin=-80,ymax=0,linestyle="dashed",color="black")
    sp.show(save=True,save_name="transmission")
