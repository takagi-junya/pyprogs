import numpy as np
import matplotlib.pyplot as plt

def input_signal(filename):
    with open(filename) as f:
        lines = f.readlines()
        signal = []
        for line in lines:
            signal.append(float(line))
    return np.array(signal)

def add_zero(signal,N):
    signal_add = np.zeros(N)
    for i in range(0,len(signal)):
        signal_add[i] = signal[i]
    return signal_add

def fft(signal):
    n = len(signal)
    fsignal = np.fft.fft(signal)
    fsignal = np.abs(fsignal/n*2)
    fsignal = fsignal[0:int(n/2)]
    return fsignal

def gausian_pulse(t,tau):
    alpha = 16.0/tau/tau
    incident = np.exp(-alpha*(t-tau)**2)
    return incident

def dgausian_pulse(t,tau):
    alpha = 16.0/tau/tau
    incident = (t-tau)/tau*np.exp(-alpha*(t-tau)**2)
    return incident
    
def plot_signal(t,signal,title="signal",\
                xlabel="time",ylabel="amplitude",\
                xmin=0,xmax=0,ymin=0,ymax=0,\
                label="signal"):
    plt.plot(t,signal,label=label)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel,rotation=90)
    if(xmin!=0 or xmax!=0):
        plt.xlim([xmin,xmax])
    if(ymin!=0 or ymax!=0):
        plt.ylim([ymin,ymax])
    
def plot_fsignal(freq,fsignal,\
                xmin=0,xmax=0,ymin=0,ymax=0,\
                color="blue",marker="",\
                xlabel="frequency",ylabel="amplitude",\
                title="fourier transformed signal",\
                label="fsignal"):
    plt.plot(freq,fsignal,color=color,marker=marker,label=label)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel,rotation=90)
    if(xmin!=0 or xmax!=0):
        plt.xlim([xmin,xmax])
    if(ymin!=0 or ymax!=0):
        plt.ylim([ymin,ymax])
    
def plot_sig_fsig(t,signal,freq,fsignal):
    fig,ax = plt.subplots(1,2,figsize=(12,4))
    ax[0].plot(t,signal)
    ax[0].grid()
    ax[1].plot(freq,fsignal)
    ax[1].grid()

def show(draw=True,save=False,save_name="signal"):
    plt.grid()
    plt.legend()
    if(save):
        plt.savefig("fig/"+save_name+".png")
    if(draw):
        plt.show()
    plt.close()
    
def rcs(finc,fscat,r=1.0,geo_area=1.0):
    return 10*np.log10(np.abs(fscat)**2/np.abs(finc)**2*r/geo_area)

def test(text,*args):
    print(text)
    print(args)