#------------------------------------------
# Plot 3D data in 1D
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import h5py
import os
#import gc
import sys
import matplotlib.patches as patches
import shutil
import make_mp4
import f90nml

params = f90nml.read("param.inp")
amp = params["scatt"]["amp"]
nx = params["space"]["nxx"]
ny = nx
#stride = params["output"]["stride"]
args = sys.argv
#set a component to draw
component = args[1]
#select drawing area 'all','center'
area = args[2]
#select type of data 'png','gif'
menu = args[3]
#plot data per 't'
t = int(args[4]) 

keys=[]
ims = []
fig,ax=plt.subplots()
cwd = os.getcwd()
os.chdir(cwd)
z0   =  376.73031

if(menu=="gif"):
  dir = cwd+"/gif"
  if not os.path.exists(dir):
      os.makedirs(dir)
    
if(menu=="png" or menu=="mp4"):
  dir = cwd+"/1dfig_"+component
  filename = "1dfig_"+component
  if not os.path.exists(dir):
      os.makedirs(dir)
  
def plot1D(data,vmin,vmax,flag,area):
  if(flag==0):
    img = ax.plot(data,color="red")
    ax.set_ylim(vmin,vmax)
  else:
    img = ax.plot(data,color="red")
  if(area==1):
    ax.set_xlim(px-nx/5,px+nx/5)
  return img

k=0
with h5py.File(component+".h5",mode="r") as f:
  for i in f[component].keys():
    keys.append(i)
  sh = f[component]["0000"].shape
  #py = int(sh[1]/2)

  for k in range(0,len(keys)):
    print(keys[k])
    if(k%t!=0):
      continue;
    lines = []
    if(menu=="png" or menu=="mp4"):
      fig,ax = plt.subplots()

    #data = np.array(f[component][str(keys[k])][py,:])
    data = np.array(f[component][str(keys[k])])
    if(component in ["ex","ey","ez"]):
      if(area=="normal"):
        #im = plot1D(data=data,vmin=-amp*0.5,vmax=amp*0.5,flag=0,area=0)
        im = plot1D(data=data,vmin=-0.001,vmax=0.001,flag=0,area=0)
      else:
        im = plot1D(data=data,vmin=-amp*0.5,vmax=amp*0.5,flag=0,area=1)
    elif(component in ["hx","hy","hz"]):
      if(area=="normal"):
        im = plot1D(data=data,vmin=-amp*0.5/z0,vmax=amp*0.5/z0,flag=0,area=0)
      else:
        im = plot1D(data=data,vmin=-amp*0.5/z0,vmax=amp*0.5/z0,flag=0,area=1)
    else:
      if(area=="normal"):
        im = plot1D(data=data,vmin=-amp*0.01,vmax=amp*0.01,flag=0,area=0)
      else:
        im = plot1D(data=data,vmin=-amp*0.01,vmax=amp*0.01,flag=0,area=1)

    plt.title(component+" on y=1280  "+str(keys[k]).zfill(4))
    lines.extend(im)
    ims.append(lines)
    plt.grid()
    if(menu=="png" or menu=="mp4"):
      plt.savefig(filename+"/"+str(keys[k]).zfill(4)+".png")
      plt.close()
  if(menu=="gif"):
    os.chdir("./gif/")
    ani = animation.ArtistAnimation(fig,ims,interval=250,blit=False)
    ani.save(component+"3d.gif",writer="imagemagick")
  if(menu=="mp4"):
    print(filename)
    os.chdir(cwd)
    os.chdir(filename)
    make_mp4.mp4(component+"_"+area)

