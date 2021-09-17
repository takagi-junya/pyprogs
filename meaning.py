import numpy as np
def mean(data,name):
    sh = data.shape
    nx = sh[1]
    ny = sh[0]
    new_data = np.zeros(sh)
    if(name in ["ex","jx","hy"]):
        for i in range(1,nx):
            for j in range(0,ny):
                new_data[i,j] = 0.5*(data[i,j]+data[i-1,j])
    elif(name in ["ey","jy","hx"]):
        for i in range(0,nx):
            for j in range(1,ny):
                new_data[i,j] = 0.5*(data[i,j]+data[i,j-1])
    elif(name in ["hz"]):
        for i in range(1,nx):
            for j in range(1,ny):
                new_data[i,j] = 0.25*(data[i,j]+data[i-1,j]+data[i,j+1],data[i,j-1])
    else:
        new_data = data
                
    return new_data