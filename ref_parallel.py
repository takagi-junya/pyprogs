import numpy as np
import matplotlib.pyplot as plt
import os
import f90nml
import analytical_solution as asl
import signal_processing as sp
import directory_edit as de
import phys_cons as pc

de.make_dir("fig")
params = f90nml.read("param.inp")
mode = int(params["scatt"]["mode"])
nx   =  params["space"]["nxx"]
dx   =  params["space"]["dx"]
io   =  params["output"]["io"]
dy   =  dx
c    =  pc.c
z0   =  pc.z0
dt   =  params["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
t0   =  params["scatt"]["tau0"]
lam  =  params["scatt"]["lambda"]*dx
prad    =  params["plasma"]["prad"]
wp = params["plasma"]["wp"]
wc = params["plasma"]["wc"][0]
wh = np.sqrt(wp**2+wc**2)
fh = wh/np.pi/2
fp = wp/2/np.pi
nu = params["plasma"]["nu"]
fc = wp/2/np.pi
pw = params["wave"]["pw"]*dx
org = params["wave"]["orgs"][0]*dx
d = 2*prad*dx
f1   =  c/(10*dx)
alpha=  16/t0/t0
amp  =  params["scatt"]["amp"]
N = params["time"]["nstep"]
N = N + 1

x = np.arange(0,nx*dx,dx)
t = np.arange(0,N*dt,dt)
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
k0  = 2*np.pi/lam

iex = sp.input_signal("exin.txt")

iey = sp.input_signal("eyin.txt")

iez = sp.input_signal("ezin.txt")

iex = sp.add_zero(iex,N)
fiex = sp.fft(iex)

iey = sp.add_zero(iey,N)
fiey = sp.fft(iey)

iez = sp.add_zero(iez,N)
fiez = sp.fft(iez)

sp.plot_signal(t,iey,xmin=dt*1000,xmax=dt*2000,color="blue",label="ey")
sp.plot_signal(t,iez,title="RCP incident wave",label="ez",color="red")
sp.show(draw=False,save=False,save_name="RCP")
sp.plot_fsignal(freq,fiey,title="RCP fourier-transformed incident wave",xmin=1e6,xmax=100e9)
sp.show(draw=False,save=False,save_name="RCPf")

rex = sp.input_signal("ex_ref.txt")
frex = sp.fft(rex)

rey = sp.input_signal("ey_ref.txt")
frey = sp.fft(rey)

rez = sp.input_signal("ez_ref.txt")
frez = sp.fft(rez)

tex = sp.input_signal("ex_trs.txt")
ftex = sp.fft(tex)

tey = sp.input_signal("ey_trs.txt")
ftey = sp.fft(tey)

tez = sp.input_signal("ez_trs.txt")
ftez = sp.fft(tez)

sp.plot_signal(t,rey,xmax=1.0e-9,title="scatterd wave")
sp.show(draw=False,save=False,save_name="scattered_wave")
sp.plot_fsignal(freq,frey,xmin=1e6,xmax=100e9,title="fscatterd_wave")
sp.show(draw=False,save=False,save_name="fscatterd_wave")

sp.plot_signal(t,tey,xmax=1.0e-9,title="transmitt wave")
sp.show(draw=False,save=False,save_name="transmit_wave")
sp.plot_fsignal(freq,ftey,xmin=1e6,xmax=100e9,title="ftransmit_wave")
sp.show(draw=False,save=False,save_name="ftransmit_wave")

fie = fiex**2+fiey**2+fiez**2
fre = frex**2+frey**2+frez**2
fte = ftex**2+ftey**2+ftez**2

'''
rcof = sp.cof(fie,fre)
tcof = sp.cof(fie,fte)
afreq,ref,trs = asl.slab_mag_parallel(wp,3e11,nu,d)
sp.plot_fsignal(freq,rcof,color="red",marker="o",xmin=1e9,xmax=100e9,ymin=0,ymax=1,label="simulation")
sp.plot_fsignal(afreq,ref,color="blue",label="analytical")
sp.show(draw=True,save=False,save_name="reflection coefficient")

sp.plot_fsignal(freq,tcof,color="red",marker="o",xmin=1e9,xmax=100e9,ymin=0,ymax=1,label="simulation")
sp.plot_fsignal(afreq,trs,color="blue",label="analytical")
sp.show(draw=True,save=False,save_name="reflection coefficient")
'''

rcs = sp.rcs(fie,fre)
tcs = sp.rcs(fie,fte)
fae = 1-sp.cof(fie,fre)-sp.cof(fie,fte)

if(mode==2):
    afreq,ref,trs = asl.slab_mag_parallelR(wp,wc,nu,d)
    afreq,ref2,trs2 = asl.slab(wp,nu,d)
    ymin = -50
    ymax = 0
    sp.plot_fsignal(freq,rcs,color="red",marker="o",xmin=1e9,xmax=100e9,ymin=ymin,ymax=ymax,label="simulation")
    sp.plot_fsignal(afreq,ref,color="blue",label="analytical(magnetized RCP)")
    sp.plot_fsignal(afreq,ref2,color="blue",title="magnetic field",ylabel="Reflection Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fp,ymin,"wp",color="black",fontsize=10)
    plt.vlines(fc,ymin=ymin,ymax=ymax,linestyle="dashed",color="green")
    plt.text(fc,ymin,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fh,ymin,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="reflection")

    sp.plot_fsignal(freq,tcs,xmin=1e9,xmax=100e9,ymin=ymin,ymax=ymax,marker="o",color="red",ylabel="Transmission Coefficient Magnitude(dB)",label="simulation")
    sp.plot_fsignal(afreq,trs,color="blue",label="analytical(magnetized RCP)")
    sp.plot_fsignal(afreq,trs2,color="blue",title="parallel magnetic field",ylabel="Transmission Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fp,ymin,"wp",color="black",fontsize=10)
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="green")
    plt.text(fc,ymin,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fh,ymin,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="transmission")

    sp.plot_fsignal(freq,fae,xmin=1e9,xmax=100e9,ymin=0,ymax=1,title="parallel magnetic filed",marker="o",color="red",ylabel="Transmission Coefficient Magnitude(dB)",label="simulation")
    #sp.plot_fsignal(afreq,trs,color="blue",label="analytical(magnetized)")
    #sp.plot_fsignal(afreq,trs2,color="blue",ylabel="Transmission Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=0,ymax=1,linestyle="dashed",color="black")
    plt.text(fp,0,"wp",color="black",fontsize=10)
    plt.vlines(fp,ymin=0,ymax=1,linestyle="dashed",color="green")
    plt.text(fc,0,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=0,ymax=1,linestyle="dashed",color="black")
    plt.text(fh,0,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="absorption")
else:
    afreq,ref,trs = asl.slab_mag_parallelL(wp,3e11,nu,d)
    afreq,ref2,trs2 = asl.slab(wp,nu,d)
    ymin = -50 
    ymax = 0
    sp.plot_fsignal(freq,rcs,color="red",marker="o",xmin=1e9,xmax=100e9,ymin=ymin,ymax=ymax,label="simulation")
    sp.plot_fsignal(afreq,ref,color="blue",label="analytical(magnetized LCP)")
    sp.plot_fsignal(afreq,ref2,color="blue",title=" magnetic field",ylabel="Reflection Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fp,ymin,"wp",color="black",fontsize=10)
    plt.vlines(fc,ymin=ymin,ymax=ymax,linestyle="dashed",color="green")
    plt.text(fc,ymin,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fh,ymin,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="reflection")

    sp.plot_fsignal(freq,tcs,xmin=1e9,xmax=100e9,ymin=ymin,ymax=ymax,marker="o",color="red",ylabel="Transmission Coefficient Magnitude(dB)",label="simulation")
    sp.plot_fsignal(afreq,trs,color="blue",label="analytical(magnetized LCP)")
    sp.plot_fsignal(afreq,trs2,color="blue",title="parallel magnetic field",ylabel="Transmission Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fp,ymin,"wp",color="black",fontsize=10)
    plt.vlines(fp,ymin=ymin,ymax=ymax,linestyle="dashed",color="green")
    plt.text(fc,ymin,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=ymin,ymax=ymax,linestyle="dashed",color="black")
    plt.text(fh,ymin,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="transmission")

    sp.plot_fsignal(freq,fae,xmin=1e9,xmax=100e9,ymin=0,ymax=1,title="parallel magnetic field",marker="o",color="red",ylabel="Transmission Coefficient Magnitude(dB)",label="simulation")
    #sp.plot_fsignal(afreq,trs,color="blue",label="analytical(magnetized)")
    #sp.plot_fsignal(afreq,trs2,color="blue",ylabel="Transmission Coefficient Magnitude(dB)",label="analytical(unmagnetized)")
    plt.vlines(fp,ymin=0,ymax=1,linestyle="dashed",color="black")
    plt.text(fp,0,"wp",color="black",fontsize=10)
    plt.vlines(fp,ymin=0,ymax=1,linestyle="dashed",color="green")
    plt.text(fc,0,"wc",color="green",fontsize=10)
    plt.vlines(fh,ymin=0,ymax=1,linestyle="dashed",color="black")
    plt.text(fh,0,"wh",color="black",fontsize=10)
    sp.show(draw=True,save=True,save_name="absorption")