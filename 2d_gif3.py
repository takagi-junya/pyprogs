#------------------------------------------
# Plot 3D data in 2D
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm
import h5py
import os
#import gc
import sys
import matplotlib.patches as patches
import shutil
import make_mp4
import f90nml
import matplotlib.patches as patches

eps0 = 8.854187817e-12
me = 9.10938188e-31
qe = 1.602176462e-19
params = f90nml.read("param.inp")
mode = params["scatt"]["mode"]
amp = params["scatt"]["amp"]
if(mode==5):
  amp = amp*0.1
  
nx = params["space"]["nxx"]
stride = params["output"]["stride"]
ny = nx
px = int((nx+params["space"]["lpml"][0]*2)/2/stride)
py = px
drad = int(params["object"]["radius"]/stride)
prad = int(params["plasma"]["prad"])
wp = params["plasma"]["wp"]
ndm =  eps0*me/(qe*qe)*wp*wp
rng = int(prad*1.1)
args = sys.argv
#set a component to draw
component = args[1]
#select drawing area 'all','center'
area = args[2]
#select type of data 'png','gif'
menu = args[3]
#plot data per 't'
t = int(args[4]) 

if(len(args)>5):
  m = int(args[5])
else:
  m = 1
  
keys=[]
ims = []
fig,ax=plt.subplots()
cwd = os.getcwd()
os.chdir(cwd)
z0   =  376.73031

#make dirctory to save image data
if(menu=="gif"):
  dir = cwd+"/gif"
  if not os.path.exists(dir):
      os.makedirs(dir)

if(menu=="png" or menu=="mp4"):
  dir = cwd+"/fig_"+component
  filename = "fig_"+component
  if os.path.exists(dir):
      shutil.rmtree(dir)
      os.makedirs(dir)
  else:
      os.makedirs(dir)

def plot2D(data,vmin,vmax,flag,area):
    if(flag==0):
      img = ax.imshow(data,vmin=vmin,vmax=vmax,cmap="jet",origin="lower")
    elif(flag==1):
      img = ax.imshow(data,cmap="jet",origin="lower")
    elif(flag==2):
      img = ax.imshow(data,vmin=ndm*0.01,vmax=ndm,cmap="jet",norm=LogNorm(),origin="lower")
    if(area==1):
      plt.xlim(px-rng,px+rng)
      plt.ylim(py-rng,py+rng)
    return img

k=0 
with h5py.File(component+".h5",mode="r") as f:
    for i in f[component].keys():
        keys.append(i) 
    sh = f[component]["0000"].shape
    x0 = int(sh[1]/2)
    y0 = int(sh[0]/2)
    pz = 0
    for k in range(0,len(keys)):
      if(k%t!=0):
        continue;  
      print(keys[k])
      data = np.array(f[component][str(keys[k])][::m,::m])
      if(component in ["ex","ey","ez"]):
        if(area=="normal"):
          im = plot2D(data=data,vmin=-amp*0.5,vmax=amp*0.5,flag=0,area=0)
        else:
          im = plot2D(data=data,vmin=-amp*0.5,vmax=amp*0.5,flag=0,area=1)
      elif(component in ["hx","hy","hz"]):
        if(area=="normal"):
          im = plot2D(data=data,vmin=-amp*0.5/z0,vmax=amp*0.5/z0,flag=0,area=0)
        else:
          im = plot2D(data=data,vmin=-amp*0.5/z0,vmax=amp*0.5/z0,flag=0,area=1)
      elif(component == "nd"):
          im = plot2D(data=data,vmin=0,vmax=0,flag=2,area=0)
      else:
        if(area=="normal"):
          im = plot2D(data=data,vmin=-amp*0.1,vmax=amp*0.1,flag=0,area=0)
        else:
          im = plot2D(data=data,vmin=-amp*0.1,vmax=amp*0.1,flag=0,area=1)
      #c = patches.Circle(xy=(x0,y0),radius=drad,fc="0.5",ec="0.5")
      #ax.add_patch(c)
      plt.xlabel("x")
      plt.ylabel("y")
      
      plt.title(component)
      if k==0:
          plt.colorbar(im)
      if(menu=="png" or menu =="mp4"):
        plt.savefig(filename+"/"+str(k).zfill(4)+".png")
      ims.append([im])
    plt.close()

if(menu=="gif"):
  os.chdir("./gif/")
  ani = animation.ArtistAnimation(fig,ims,interval=250,blit=True)
  ani.save(component+"_"+area+"3d.gif",writer="imagemagick")
elif(menu=="mp4"):
  os.chdir(filename+"/")
  make_mp4.mp4(component+"_"+area)
