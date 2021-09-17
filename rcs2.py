import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import solve_te
import f90nml

cwd = os.getcwd()
dir = cwd+"/fig"
os.chdir(cwd)
if not os.path.exists(dir):
    os.makedirs(dir)

params = f90nml.read("param.inp")

keys = []
ff = []
k = 0

dx   =  params["space"]["dx"]
dy   =  dx
c    =  2.9979246e8
z0   =  376.73031
dt   =  params["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
t0   =  params["scatt"]["tau0"]
a    =  params["object"]["radius"]
a    =  a*dx
f1   =  c/(10*dx)
alpha=  16/t0/t0
amp  =  params["scatt"]["amp"]
wp   =  params["plasma"]["wp"]

N = params["time"]["nstep"]
N = N 
t = np.arange(0,N*dt,dt)
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
k0 = 2*np.pi*freq/c
print("fmax:,{:.3g}".format(f1))

ez = np.array(amp*np.exp(-alpha*(t-t0)**2))
ez = np.array(amp*(t-t0)/t0*np.exp(-alpha*(t-t0)**2))

with h5py.File("dz.h5",mode="r") as f1:
  dz = f1["dz"]["dz"][0:N]
with h5py.File("dphi.h5",mode="r")   as f2:
  dphi = f2["dphi"]["dphi"][0:N]

fez = np.fft.fft(ez)
fez = np.abs(fez/N*2)
fez = fez[0:int(N/2)]

fdz = np.fft.fft(dz)
fdz = np.abs(fdz/N*2)
fdz = fdz[0:int(N/2)]

fdphi = np.fft.fft(dphi)
fdphi = np.abs(fdphi/N*2)
fdphi = fdphi[0:int(N/2)]

fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,dz)
ax[0].set_title("Dz")
ax[0].grid()
ax[1].plot(freq,fdz[0:N])
#ax[1].set_xlim([0,12])
ax[1].set_title("fft_Dz")
ax[1].grid()
ax[1].set_xlim([0,20e9])
plt.savefig("fig/fdz.png")
plt.show()
plt.close()

fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,dphi)
ax[0].set_title("Dphi")
ax[0].grid()
ax[1].plot(freq,fdphi[0:N])
ax[1].set_xlim([0,20e9])
ax[1].set_title("fft_Dphi")
ax[1].grid()
plt.savefig("fig/fdphi.png")
plt.show()
plt.close()

fdz = np.abs(fdz)**2
fdphi = np.abs(fdphi)**2
fez  = np.abs(fez)**2

ka2 = np.linspace(0.1,12,1000)
sigma = k0*(fdphi+fdz)/fez/4
rcs = 10*np.log10(sigma/np.pi/a)
rcs2= solve_te.TE(ka2,a)

freq2 = ka2*c/(2*np.pi*a)
fig,ax = plt.subplots()
plt.plot(freq,rcs,"-",color="b",label="simliation")
plt.plot(freq2,rcs2,color="r",label="theory")
plt.ylabel("RCS (dB)")
plt.xlim([0,20e9])
#plt.ylim([-20,5])
#plt.xlim([0,12])
plt.legend()
plt.grid()
plt.savefig("fig/rcs2.png")
plt.show()
plt.close()