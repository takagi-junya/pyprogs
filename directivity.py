from typing import TYPE_CHECKING
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
import f90nml

args = sys.argv

menu = args[1]
print(menu)
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
nu   =  params["plasma"]["nu"]
prad =  params["plasma"]["prad"]*dx
cfreq = freq
lm= c/cfreq
comega=2*np.pi*cfreq
tau0 = params["scatt"]["tau0"]
ntau0=int(tau0/dt)   
radi0 = 1.74532925e-2
phi = radi0*phi0
r0x = np.cos(phi)
r0y = np.sin(phi)
vbc = c
wpw = wp/comega
nuw  = nu/comega
cr = 2*(prad-a)/lm

cwd = os.getcwd()
dir = cwd+"/fig"
os.chdir(cwd)
if not os.path.exists(dir):
    os.makedirs(dir)

cwdd = cwd.split("/")
dirl = float(cwdd[len(cwdd)-1])
with h5py.File("D.h5",mode="r") as f:
    D = np.array(f["D"]["D"])
with h5py.File("/LARGE0/gr20001/b36288/pythons/fdtd/Db.h5",mode="r") as f:
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
rcs = 10*np.log10(sigma/np.pi/a)

yb = Db/np.abs(E0)**2
sigmab = k*yb/4
rcsb = 10*np.log10(sigmab/np.pi/a)

#rcs = rcsb - rcs
import out_sl
wpw = out_sl.round(wpw)
nuw = out_sl.round(nuw)
cr = out_sl.round(cr)
if(menu=="nu"):
  xl_name = "nu.xlsx"
  tag = "nu/w"
if(menu=="wp" or menu=="wpg"):
  xl_name = "wp.xlsx"
  tag = "wp/w"
if(menu=="rp"):
  xl_name = "rp.xlsx"
  tag = "rp"

xl_path = out_sl.parentDir()+"/"+xl_name
print("wp/w:",wpw)
print("nu/w:",nuw)
print("2(rp-rc)/lm:",cr)
book = out_sl.init(xl_path,xl_name,tag)
sheet = book["directivity"]
if(menu=="nu"):
  out_sl.add_data(sheet,out_sl.find(sheet,nuw),nuw,rcs)
if(menu=="wp"):
  out_sl.add_data(sheet,out_sl.find(sheet,wpw),wpw,rcs)
if(menu=="wpg"):
  out_sl.add_data(sheet,out_sl.find(sheet,dirl),dirl,rcs)
if(menu=="rp"):
  out_sl.add_data(sheet,out_sl.find(sheet,cr),cr,rcs)

book.save(xl_path)
'''
ax = plt.subplot(111,projection="polar")
ax.plot(x,rcs,"-",label="with plasma",linewidth=2,color="red")
if(menu=="nu"):
  ax.set_title("RCS(dB) "+"ν/$\omega$="+str(nuw)+" "+"2$(r_p-r_c)$/$\lambda$=0.1",pad=10)
if(menu=="wp"):
  ax.set_title("RCS(dB) "+"$\omega_p$/$\omega$="+str(wpw)+" "+"2$(r_p-r_c)$/$\lambda$=0.1",pad=10)
if(menu=="rp"):
  ax.set_title("RCS(dB) "+"$\omega_p$/$\omega$="+str(wpw)+" "+"2$(r_p-r_c)$/$\lambda$="+str(cr),pad=10)
ax.set_xlabel("observation angle")
plt.legend()
os.chdir(cwd)
plt.savefig("fig/dir.png")
#plt.show()
plt.close()
'''