import sys
import os 
import time
import subprocess
import numpy as np 
import out_sl

cwd = os.getcwd()
args = sys.argv 
menu = args[1]

files = os.listdir(cwd)
files_dir = [direc for direc in files if os.path.isdir(cwd+"/"+direc)]

pathes = []
for dirc in files_dir:
    pathes.append(cwd+"/"+dirc)

i = 0
for path in pathes:
    os.chdir(path)
    print(files_dir[i])
    subprocess.run(["/LARGE0/gr20001/b36288/python3","/LARGE0/gr20001/b36288/pythons/fdtd/directivity.py",menu])
    i = i + 1
    time.sleep(5)
    
os.chdir(cwd)
xl_name = menu+".xlsx"
if(menu=="wp" or menu=="nu"):
    tag = menu+"/w"
elif(menu=="wpg"):
    tag = "wp/w"
    menu = "wp"
elif(menu=="rp"):
    tag = menu
out_sl.sort(xl_name,tag)