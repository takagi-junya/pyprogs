import f90nml
import sys
import os 
import shutil
import numpy as np
import subprocess
import integral2d
import math
from scipy.special import *
    
cwd = os.getcwd()
args = sys.argv

root = args[1]
dest = args[2]
opt  = args[3]
coef = float(dest)
fdtd_path = "/home/b/b36288/large0/drude/fdtd"

root_path = cwd+"/"+root
dest_path = cwd+"/"+dest

if not os.path.exists(dest_path):
    os.mkdir(dest_path)
    
if os.path.exists(dest_path+"/fdtd"):
    os.remove(dest_path+"/fdtd")
    
os.chdir(dest_path)
shutil.copyfile(root_path+"/job.sh",dest_path+"/job.sh")
os.symlink(fdtd_path,dest_path+"/fdtd")

with open(root_path+"/param.inp","r") as nml_file:
    params = f90nml.read(nml_file)
    
dest_nml = params
with open(dest_path+"/param.inp","w") as nml_file:
    f90nml.write(dest_nml,dest_path+"/param.inp",force=True)
    
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
freq = params["scatt"]["freq"]
phi0 =  params["scatt"]["phi0"]
wp   =  params["plasma"]["wp"]
nu   =  params["plasma"]["nu"]
prad = params["plasma"]["prad"]*dx
cfreq = freq
lm = c/cfreq
comega=2*np.pi*cfreq
k0a = comega*a/c
tau0 = params["scatt"]["tau0"]
ntau0=int(tau0/dt)   
radi0 = 1.74532925e-2
phi = radi0*phi0
r0x = np.cos(phi)
r0y = np.sin(phi)
vbc = c
wpw = wp/comega
nuw  = nu/comega
    
def wpg():
    px = int(nx/2)
    py = px
    ds = dx*dx
    rad = int(a/dx)
    pprad = int(prad/dx)
    def radi(i,j):
        return np.sqrt(((i-px))**2+((j-py))**2)

    count = 0
    nancount = 0
    n = np.zeros([nx,nx])
    lower = int(px-pprad*1.05)
    upper = int(px+pprad*1.05)
    for i in range(lower,upper):
        for j in range(lower,upper):
            if(radi(i,j)>=rad and radi(i,j)<pprad):
                n[i,j] = wp*np.sqrt(jn(0,2.40*(radi(i,j)-rad)/(pprad-rad)))
                if(math.isnan(n[i,j])):
                    n[i,j] = 0.0
                    nancount = nancount + 1
                count = count + 1

    dn = integral2d.integ2d(n,dx)
    nds = ds*count
    print("wp/w={:3f}".format(dn/nds/comega))
    
if(opt=="wp"):
    wp = comega*coef
    dest_nml["plasma"]["wp"] = wp
    
if(opt=="nu"):
    nu = comega*coef
    dest_nml["plasma"]["nu"] = nu
    
if(opt=="wpg"):
    #rad50#wp = comega*coef*1.38555552
    #wp = comega*coef*1.42467882
    #rad100#wp = comega*coef*1.39210827
    wp = comega*coef*1.42779182
    dest_nml["plasma"]["wp"] = wp
    dest_nml["plasma"]["pls"] = 2

if(opt=="rp"):
    prad = int((lm*coef/2+a)/dx)
    dest_nml["plasma"]["prad"] = prad
    prad = prad*dx
    
with open(dest_path+"/param.inp","w") as nml_file:
    f90nml.write(dest_nml,dest_path+"/param.inp",force=True)


print("freq:{:3e}".format(cfreq))
print("wp:{:3e}".format(wp))
print("nu:{:3e}".format(nu))
print("nx:",nx)
print("rad:",a)
print("prad:",prad)
print("lambda:",lm)
print("omega:{:3e}".format(comega))
print("ka:{:.3f}".format(k0a))
if(opt=="wpg"):
    wpg()
else:
    print("wp/w:{:.3f}".format(wp/comega))
print("nu/w:{:.3f}".format(nu/comega))
print("2(rp-rc)/lm:{:.3f}".format(2*(prad-a)/lm))

subprocess.run(["qsub","job.sh"])

'''
nml = {
    "space":{
        "nxx":3000,
        "nyy":3000,
        "dx":0.00025,
        "dy":0.00025,
        "abc":1,
        "pbc":[0,0,0],
        "lpml":[8,8,0]
    },
    "time":{
        "deltat":0.80,
        "nstep":12000
    },
    "output":{
        "out":20000,
        "ostart":0,
        "odom":[0,0,0,0],
        "stride":1,
        "comp":[0,0,0,0,0,0,0,0,0],
        "io":600,
        "jo":1008
    },
    "scatt":{
        "mode":3,
        "lx":450,
        "ly":450,
        "gamma0":-90.0,
        "thetat0":90.0,
        "phi0":180.0,
        "amp":1.0,
        "lambda":400,
        "tau0":1.0e-20
    },
    "far":{
        "isx":5,
        "isy":5,
        "theta1":90.0,
        "phi1":180.0,
    },
    "object":{
        "obj":1,
        "med":2,
        "ic":508,
        "jc":508,
        "lx2":0,
        "ly2":0,
        "epsr":1.0,
        "radius":40.0,
    },
    "feed":{
        "lfeed":0,
        "ip":0,
        "jp":0,
        "duration":0.0,
        "t0":0.0
    },
    "wave":{
        "kwave":0,
        "amps":[0.0,0.0,1.0,0.0,1.0,0.0],
        "orgs":[150.0,0.0,0.0],
        "angs":[90.0,0.0,0.0],
        "pt":5.0,
        "pw":10.0
    },
    "plasma":{
        "pls":1,
        "prad":60,
        "nu":0.0,
        "wp":3.3903e10,
    }
}
'''