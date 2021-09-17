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

nx   =  params["space"]["nxx"]
ny   =  params["space"]["nyy"]
px   =  int(nx/2)
py   =  int(ny/2)
dx   =  params["space"]["dx"]
dy   =  dx
c    =  2.9979246e8
z0   =  376.73031
dt   =  params["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
tau0 =  params["scatt"]["tau0"]
lx   =  params["scatt"]["lx"]
a    =  params["object"]["radius"]*dx
cfreq = params["scatt"]["freq"]
comega= 2*np.pi*cfreq
alpha=  16/tau0/tau0
amp  =  params["scatt"]["amp"]
wp   =  params["plasma"]["wp"]

N = params["time"]["nstep"]
N = N
t = np.arange(0,N*dt,dt)
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
k0 = 2*np.pi*freq/c
for index in range(0,len(freq)):
  if(freq[index]>3.0e9):
    break
  index = index+ 1

print("freq:{:.3e},".format(freq[index]))
tau = t
tt  = tau-tau0
einc = np.array(amp*(tt)/tau0*np.exp(-alpha*tt*tt))
plt.plot(einc)
plt.grid()
plt.show()

feinc = np.fft.fft(einc)
feinc = np.abs(feinc/N*2)
feinc = feinc[0:int(N/2)]
plt.plot(freq,feinc)
plt.xlim([0,10e9])
plt.grid()
plt.show()

with h5py.File("dz.h5",mode="r") as f1:
  dz = f1["dz"]["dz"][0:N]
with h5py.File("dphi.h5",mode="r")   as f2:
  dphi = f2["dphi"]["dphi"][0:N]

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
ax[1].set_xlim([0,10e9])
plt.savefig("fig/fdz.png")
plt.show()
plt.close()

fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,dphi)
ax[0].set_title("Dphi")
ax[0].grid()
ax[1].plot(freq,fdphi[0:N])
ax[1].set_xlim([0,10e9])
ax[1].set_title("fft_Dphi")
ax[1].grid()
plt.savefig("fig/fdphi.png")
plt.show()
plt.close()

fdz = np.abs(fdz[0:N])**2
fdphi = np.abs(fdphi[0:N])**2
feinc  = np.abs(feinc[0:N])**2

#ka2 = np.linspace(0.1,12,1000)
range = np.where(freq<5.5e9)
freq = freq[range]
ka2 = 2*np.pi*freq*a/c
sigma = k0[0:N]*(fdphi+fdz)/feinc/4
rcs = 10*np.log10(sigma[range]/np.pi/a)
rcs2= solve_te.TE(ka2,a)

freq2 = ka2*c/(2*np.pi*a)
fig,ax = plt.subplots()
plt.plot(freq,rcs,"-",color="b",label="simliation")
plt.plot(freq,rcs2,color="r",label="theory")
plt.ylabel("RCS (dB)")
plt.xlim([0,5e9])
plt.legend()
plt.grid()
plt.show()
plt.close()

print("vacum:{:.3e},plasma:{:.3e}".format(rcs2[index],rcs[index]))
print("diff:{:.3f}".format(rcs[index]-rcs2[index]))

'''
'''