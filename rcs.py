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
print(N)
t = np.arange(0,N*dt,dt)
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
k0  = 2*np.pi*freq/c
ka = k0*a
kp = wp/c*a
#ka=ka[np.where(ka<12)]

print("fmax:,{:.3g}".format(f1))
print("kamax:",2*np.pi*a/(10*dx))

ez = np.array(amp*np.exp(-alpha*(t-t0)**2))

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

'''
fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,ez)
ax[0].set_title("sig_ez")
ax[0].grid()
ax[1].plot(ka,fez[0:N])
ax[1].set_title("fft_ez")
ax[1].grid()
ax[1].set_xlim([0,12])
plt.savefig("fig/ez.png")
plt.show()
plt.close()

fig,ax = plt.subplots()
ax.plot(freq,fez)
ax.grid()
plt.show()
plt.close()
'''

fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,dz)
ax[0].set_title("Dz")
ax[0].grid()
ax[1].plot(ka,fdz[0:N])
ax[1].set_xlim([0,12])
ax[1].set_title("fft_Dz")
ax[1].grid()
ax[1].set_xlim([0,12])
plt.savefig("fig/fdz.png")
plt.show()
plt.close()

fig,ax = plt.subplots(2,1,figsize=(5,8))
ax[0].plot(t,dphi)
ax[0].set_title("Dphi")
ax[0].grid()
ax[1].plot(ka,fdphi[0:N])
ax[1].set_xlim([0,12])
ax[1].set_title("fft_Dphi")
ax[1].grid()
plt.savefig("fig/fdphi.png")
plt.show()
plt.close()

fdz = np.abs(fdz[0:N])**2
fdphi = np.abs(fdphi[0:N])**2
fez  = np.abs(fez[0:N])**2

n1 = np.argmax(fdphi)
print(ka[n1])

ka2 = np.linspace(0.20,12,2000)
sigma = k0[0:N]*(fdphi+fdz)/fez/4
rcs = 10*np.log10(sigma/np.pi/a)
rcs2= solve_te.TE(ka2,a)

fig,ax = plt.subplots()
plt.plot(ka,rcs,"-",color="b",label="simliation")
plt.plot(ka2,rcs2,color="r",label="theory")
plt.xlim([0,12])
plt.ylim([-10,5])
plt.vlines(kp,-10,5,color="red",linestyle="dashed")
plt.text(kp,3,"$\omega_p$",fontsize=15)
plt.xlabel("ka="+"2$\pi$a/$\lambda$")
plt.ylabel("RCS (dB)")
plt.title("$\omega_p=14.0 GHz$")
plt.legend()
plt.grid()
plt.savefig("fig/rcs.png")
plt.show()
plt.close()

'''
plt.plot(t,inc)
plt.grid()
plt.show()

plt.plot(t,scat)
plt.grid()
plt.show()
'''

'''
finc = np.fft.fft(inc)
finc = np.abs(finc/N*2)
finc = finc[0:int(N/2)]

fscat = np.fft.fft(scat)
fscat = np.abs(fscat/N*2)
fscat = fscat[0:int(N/2)]

r = 700+8
r = r*dx
sigma3 = 2*np.pi*(r-a/2)*np.abs(fscat)**2/np.abs(finc)**2
rcs3 = 10*np.log10(sigma3/np.pi/a)
plt.plot(ka,rcs3,label="near-field",color="black")
plt.plot(ka,rcs,label="far-field",color="red")
plt.plot(ka2,rcs2,label="theory",color="blue")
plt.xlim([0,5])
plt.ylim([-3,3])
plt.grid()
plt.legend()
plt.show()
'''