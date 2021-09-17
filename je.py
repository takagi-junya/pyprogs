import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import animation
import f90nml
import os
import h5py
import meaning
import integral2d

param = f90nml.read("param.inp")
cwd = os.getcwd()
dx = param["space"]["dx"]
def meanh5(tag,k):
    with h5py.File(tag+".h5",mode="r") as f:
        data = np.array(f[tag][str(k).zfill(4)])
        meaning.mean(data,tag)
    return data

fig,ax = plt.subplots()
with h5py.File("ex.h5",mode="r") as f:
    keys = list(f["ex"].keys())
i = 0
int_je = 0.0
int_jel = []
ims = []
for k in (keys):
    print(k)
    ex = meanh5("ex",i)
    ey = meanh5("ey",i)
    ez = meanh5("ez",i)
    jx = meanh5("jx",i)
    jy = meanh5("jy",i)
    jz = meanh5("jz",i)
    
    je = ex*jx+ey*jy+ez*jz
    sje = integral2d.integ2d(je,dx)
    int_je = int_je + sje
    int_jel.append(sje)
    
    im = plt.imshow(je,vmin=-0.1,vmax=0.1,cmap="jet",origin='lower')
    if(i==0):
        plt.title("JE")
        plt.colorbar(im)
    ims.append([im])
    i = i + 1
os.chdir("./gif/")
ani = animation.ArtistAnimation(fig,ims,interval=250,blit=True)
ani.save("je.gif",writer="imagemagick")
plt.close()

'''
os.chdir(cwd)
fig,ax = plt.subplots(figsize=(12,4))
print(int_je)
plt.plot(int_jel,label="JE")
plt.title("JE")
plt.legend(loc="upper left")
plt.grid()
plt.savefig("fig/je.png")
plt.show()
plt.close()
'''
