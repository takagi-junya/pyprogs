import numpy as np 
import matplotlib.pyplot as plt
import f90nml
import integral2d
from scipy.special import *
import math


c = 2.99792458e8
param = f90nml.read("param.inp")
dx = param["space"]["dx"]
lpml = param["space"]["lpml"][0]
dt   =  param["time"]["deltat"]/(c*np.sqrt(1.0/(dx*dx)+1.0/(dx*dx)))
nx = param["space"]["nxx"]
rad = param["object"]["radius"]
prad = param["plasma"]["prad"]
nu =param["plasma"]["nu"]
wp = param["plasma"]["wp"]
#lam = param["scatt"]["lambda"]
freq = param["scatt"]["freq"]
pls = int(param["plasma"]["pls"])
a = rad*dx

#lm = lam*dx
lm = c/freq
lam = lm/dx
#freq = c/lm
omega = 2*np.pi*freq
k0a = omega*a/c
r = prad-rad
ds = dx*dx
print("CFL:{:3f}".format(np.sqrt(1-(wp*dt/2)**2)))
print("freq:{:3e}".format(freq))
print("wp:{:3e}".format(wp))
print("nx:",nx)
print("rad:",rad)
print("prad:",prad)
print("lambda:",lm)
print("omega:{:3e}".format(omega))
print("ka:{:.3f}".format(k0a))

if(pls>1):
    print("wp/w:{:.3f}".format(wp/omega))

'''
elif (pls==2):
    px = int(nx/2)
    py = px
    def radi(i,j):
        return np.sqrt(((i-px))**2+((j-py))**2)

    count = 0
    nancount = 0
    wp = wp * omega
    n = np.zeros([nx,nx])
    lower = int(px-prad*1.2)
    upper = int(px+prad*1.2)
    for i in range(lower,upper):
        for j in range(lower,upper):
            if(radi(i,j)>=rad and radi(i,j)<prad):
                n[i,j] = wp*np.sqrt(jn(0,2.41*(radi(i,j)-rad)/(prad-rad)))
                if(math.isnan(n[i,j])):
                    n[i,j] = 0.0
                    nancount = nancount + 1
                count = count + 1

    dn = integral2d.integ2d(n,dx)
    ds = ds*count
    print("wp={:.3e}".format(dn/ds))
    print("wp/w={:3f}".format(dn/ds/omega))
    plt.plot(n[:,px],"o")
    plt.xlim([px-prad*1.1,px+prad*1.1])
    plt.grid()
    #plt.show()
    '''
print("nu/w:{:.3f}".format(nu/omega))
print("2(rp-rc)/lm:{:.3f}".format(2*(prad-rad)/lam))

data = []
with open("wp.txt",mode="r") as f:
    lines = f.readlines()
    for line in lines:
        data.append(float(line))
del data[0:lpml]
del data[-lpml:]
plt.plot(data)
plt.xlim([1200,1800])
plt.grid()
plt.show()