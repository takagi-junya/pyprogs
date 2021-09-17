import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import f90nml
import h5py
import phys_cons as pc
import os 
import sys
import directory_edit as DE

c = pc.c
z0 = pc.z0
class Visualizer:
    def __init__(self):
        self.params = f90nml.read("param.inp")
        self.nx = self.params["space"]["nxx"]
        self.ny = self.params["space"]["nyy"]
        self.stride=1# self.params["output"]["stride"]
        self.cx = int(self.nx/2/self.stride)
        self.cy = int(self.ny/2/self.stride)
        self.rad = int(self.params["object"]["radius"]/self.stride)
        self.prad= int(self.params["plasma"]["prad"]/self.stride)
        self.amp = self.params["scatt"]["amp"]
        self.fig,self.ax =plt.subplots(figsize=(8,4))
        
    def readh5(self,filename,tag):
        with h5py.File(filename,mode="r") as f:
            keys = list(f[tag])
            sh = f[tag][keys[0]].shape
        return keys,sh
        
    def plot_1d(self,data1d,color="blue",xmin=0,xmax=0,ymin=0,ymax=0,\
                title="",xlabel="",ylabel="",\
                ):
        img = self.ax.plot(data1d,color=color)
        if(xmax==0):
            xmax = int(data1d.shape[0])
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_xlim([xmin,xmax])
        self.ax.set_ylim([ymin,ymax])
        return img
    
    def plot_2d(self,data2d,vmin,vmax,xmin=0,xmax=0,ymin=0,ymax=0,\
                title="",xlabel="",ylabel="",\
                ):
        img = self.ax.imshow(data2d,vmin=vmin,vmax=vmax,cmap="jet",origin="lower")
        if(xmax==0):
            xmax = data2d.shape[1]
        if(ymax==0):
            ymax = data2d.shape[0]
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_xlim([xmin,xmax])
        self.ax.set_ylim([ymin,ymax])
        return img
    
    def movie_1d(self,filename,tag,cof=0.5,start=0,interval=1,figsave=False,savename="fig",tsave="gif"):
        ims = []
        DE.make_dir(tag+"fig")
        DE.make_dir("movie")
        keys,sh = self.readh5(filename,tag)
        
        if(tag in ["ex","ey","ez"]):
            vmin = -self.amp*cof
            vmax = -vmin
        elif(tag in ["hx","hy","hz"]):
            vmin = -self.amp*cof/pc.z0
            vmax = -vmin
        else:
            vmin = -self.amp*cof
            vmax = -vmin
            
        if(tsave=="mp4"):
            figsave = True
            
        with h5py.File(filename,mode="r") as f:
            for i in range(start,len(keys),interval):
                print(keys[i])
                data = f[tag][keys[i]][int(sh[0]/2),:]
                im = self.plot_1d(data,ymin=vmin,ymax=vmax,title=tag,ylabel=tag)
                ims.append(im)
                plt.grid()
                if(figsave):
                    self.save_png(tag+"fig/"+savename+str(keys[i]).zfill(4))
            if(tsave=="gif"):
                ani = animation.ArtistAnimation(self.fig,ims,interval=250,blit=False)
                ani.save("movie/"+tag+"1d.gif",writer="imagemagick")
            elif(tsave=="mp4"):
                self.save_mp4(tag)
                
    def movie_2d(self,filename,tag,cof=0.5,start=0,interval=1,figsave=False,savename="fig",tsave="gif"):
        ims = []
        DE.make_dir(tag+"2dfig")
        DE.make_dir("movie")
        keys,sh = self.readh5(filename,tag)
        
        if(tag in ["ex","ey","ez"]):
            vmin = -self.amp*cof
            vmax = -vmin
        elif(tag in ["hx","hy","hz"]):
            vmin = -self.amp*cof/pc.z0
            vmax = -vmin
        else:
            vmin = -self.amp*cof
            vmax = vmin
            
        if(tsave=="mp4"):
            figsave = True
            
        with h5py.File(filename,mode="r") as f:
            
            for i in range(start,len(keys),interval):
                print(keys[i])
                data = f[tag][keys[i]][:,:]
                im = self.plot_2d(data,vmin=vmin,vmax=vmax,title=tag,xlabel="x",ylabel="y")
                ims.append(im)
                
                if(i==start):
                    plt.colorbar(im)
                if(figsave):
                    self.save_png(tag+"2dfig/"+savename+str(keys[i].zfill(4)))
            
            if(tsave=="gif"):
                ani = animation.ArtistAnimation(self.fig,ims,intercal=250,blit=False)
                ani.save("movie/"+tag+"2d.gif",write="imagemagick")
            elif(tsave=="mp4"):
                self.save_mp4(tag+"2d")
        
    def save_png(self,savename="fig"):
        plt.savefig(savename+".png")
        plt.close()
        self.fig,self.ax = plt.subplots()
    
    def save_mp4(self,savename):
        import cv2
        import glob
        import pathlib
        cwd = os.getcwd()
        frame_rate=3.0
        width = 640
        height= 480
        fourcc = cv2.VideoWriter_fourcc('m','p','4','v')
        video = cv2.VideoWriter(savename+".mp4",fourcc,frame_rate,(width,height))
        print("making mp4 file")
        os.chdir(savename+"fig")
        images = sorted(glob.glob("*png"))
        for i in range(len(images)):
            img = cv2.imread(images[i])
            img = cv2.resize(img,(width,height))
            video.write(img)
        video.release()
        os.chdir(cwd)
        DE.remove_dir(savename+"fig")
        
print()
dim = input("dimension 1d or 2d:")
print()
tag = input("output component:")
print()
interval = int(input("output interval:"))
print()
tsave = input("save file type:")

if(dim=="1d"):
    vs = Visualizer()
    vs.__init__()
    vs.movie_1d(filename=tag+".h5",cof=1.5,start=0,interval=interval,figsave=False,savename=tag,tag=tag,tsave=tsave)
elif(dim=="2d"):
    vs = Visualizer()
    vs.__init__()
    vs.movie_2d(filename=tag+".h5",cof=0.5,start=0,interval=interval,figsave=False,savename=tag,tag=tag,tsave=tsave)