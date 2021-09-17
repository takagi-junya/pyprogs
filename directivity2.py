import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import f90nml
params = f90nml.read("param.inp")
dx   =  params["space"]["dx"]
nx   =  params["space"]["nxx"]
lx   =  params["scatt"]["lx"]
dy   =  dx
c    =  2.9979246e8
dt   =  params["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
nstep=  params["time"]["nstep"]
a    =  params["object"]["radius"]*dx
#lm   =  params["scatt"]["lambda"]*dx
amp = params["scatt"]["amp"]
freq = params["scatt"]["freq"]
phi0 =  params["scatt"]["phi0"]
wp   =  params["plasma"]["wp"]
lm = c/freq
cfreq= freq
comega=2*np.pi*cfreq
tau0 = params["scatt"]["tau0"]
ntau0=int(tau0/dt)   
radi0 = 1.74532925e-2
phi = radi0*phi0
r0x = np.cos(phi)
r0y = np.sin(phi)
vbc = c
wpw = wp/comega
print("dx:",dx)
print("a:",a)
print("dt:",dt)
print("comega:{:.3e}".format(comega))

cwd = os.getcwd()
dir = cwd+"/fig"
os.chdir(cwd)
if not os.path.exists(dir):
    os.makedirs(dir)

with h5py.File("D.h5",mode="r") as f:
    D = np.array(f["D"]["D"])

''' 
with h5py.File("D2.h5",mode="r") as f:
    D2 = np.array(f["D"]["D"])
'''
'''
with h5py.File("/LARGE0/gr20001/b36288/pythons/fdtd/Db.h5",mode="r") as f:
    Db = np.array(f["D"]["D"])
'''

with h5py.File("Db.h5",mode="r") as f:
    Db = np.array(f["D"]["D"])
i = 0

t=0.0
N = nstep+1
inc = np.zeros([N])
ny = nx
qx = nx*dx*r0x
qy = ny*dy*r0y
dis = 0.0 
dd = qx
if(dd>dis): dis = dd
dd = qy 
if(dd>dis): dis = dd
dd = qx+qy
if(dd>dis): dis = dd
x = (lx+0)*dx
y = (int(ny/2)+0.5)*dy

for i in range(0,N):
  t = t+dt
  tau = t+(r0x*x+r0y*y-dis)/c
  tt = tau-tau0
  if(tt<2.0*np.pi/comega):
    w = 0.0
  elif(tt>=2.0*np.pi/comega and tt<=3.0*np.pi/comega):
    w = 0.5*(1.0-np.cos(comega*tt))
  else:
    w = 1.0
  inc[i] = amp*w*np.sin(comega*tt)
  
finc = np.fft.fft(inc)
finc = np.abs(finc/N*2)
finc = finc[0:int(N/2)]
freq = np.fft.fftfreq(N,dt)
freq = freq[0:int(N/2)]
plt.plot(freq,finc)
plt.grid()
plt.xlim([0,10e9])
plt.show()

nintg = 3
tintg = nintg/cfreq
ndt = int(tintg/dt)+1
if(ndt%2!=0.0):
   ndt = ndt+1
itgstr = nstep-ndt
itgend = nstep
wi = 1/3*dt/tintg

cj = 0.0+1.0j
co = 0.0+0.0j
for step in range(itgstr,itgend+1):
  t = step*dt
  cexpe = np.exp(-cj*comega*t)
  if(step==itgstr or step==itgend):
    cof = wi*cexpe
  elif(step%2==0):
    cof = 2.0*wi*cexpe
  else:
    cof = 4.0*wi*cexpe
  co = co + cof*inc[step]
E0 = np.abs(co)
print("E0:",E0)
t = np.arange(0,N)
k = comega/c 
ka=k*a

p = 72
d = 5.0

x = np.arange(0,360+d,d)
x1 = x
x = (np.pi/180)*x

y = D/np.abs(E0)**2
sigma = k*y/4
rcs = 10*np.log10(sigma/(2.4*0.1))#/np.pi/a)
rcs = 10*np.log10(sigma/lm/2)
'''
y2 = D2/np.abs(E0)**2
sigma2 = k*y2/4
rcs2 = 10*np.log10(sigma2/np.pi/a)
'''

yb = Db/np.abs(E0)**2
sigmab = k*yb/4.0
rcsb = 10*np.log10(sigmab/lm/2)#np.pi/a)
#rcsb =rcs

#rcsb = 10*np.log10(sigmab/b)
#print(rcsb)

#for i in range(0,int(len(rcs)/2),1):
#  print("{:.0f} {:.3e}".format(x1[i],rcs[i]-rcsb[i]))

'''
for i in range(0,int(len(rcs)/2),1):
  print(x1[i],rcs[i])
  print(x1[len(x1)-1-i],rcs[len(rcs)-1-i])
  print(rcs[i]-rcs[len(rcs)-1-i])
  print()
print(x1[int(len(rcs)/2)],rcs[int(len(rcs)/2)])
'''

'''
for i in range(0,int(len(rcsb)/2),1):
  print(x1[i],rcsb[i])
  print(x1[len(x1)-1-i],rcsb[len(rcsb)-1-i])
  print(rcsb[i]-rcsb[len(rcsb)-1-i])
  print()

print(x1[int(len(rcsb)/2)],rcsb[int(len(rcsb)/2)])
'''

'''
with open("dir.txt",mode="w") as f:
  f.write("kka:"+str(kka)+"\n")
  for i in range(0,len(rcs),3):
    f.write(str(x1[i])+","+str(rcs[i])+"\n")
'''

ax = plt.subplot(111,projection="polar")
ax.plot(x,rcsb,"-",label="no plasma",linewidth=2,color="blue")
ax.plot(x,rcs,"-",label="no grad",linewidth=2,color="red")
#ax.plot(x,rcs2,"-",label="with grad",linewidth=2,color="green")
#ax.set_title("RCS(dB) "+"$\omega_p$/$\omega$=5.0 "+"2$(r_p-r_c)$/$\lambda$=0.1",pad=10)
#ax.set_title("RCS(dB) "+"ν/$\omega$=2.0 "+"2$(r_p-r_c)$/$\lambda$=0.1",pad=10)
ax.set_xlabel("observation angle")
ax.legend()
#plt.savefig("fig/dir.png")
plt.show()