import numpy as np
import matplotlib.pyplot as plt

def slab_th(wp,nu,d):
    c  = 2.998e8      #光速
    freq = np.linspace(1e6,80e9,1000) #周波数
    w  = 2*np.pi*freq #角周波数
    k0 = 2*np.pi*freq/c #波数
    
    ref   = np.zeros([len(freq)])   #反射強度
    trs   = np.zeros([len(freq)])   #透過強度

    print("wp:{:.3e}".format(wp))
    print("nu:{:.3e}".format(nu))
    print("d:{:.3f}".format(d))
    
    for i in range(0,len(freq)):
        a = wp**2*w[i]
        b = w[i]*(w[i]**2+nu**2)
        alpha = 1 - a/b  #比誘電率の実部
        beta  = -wp**2*nu/(b)   #比誘電率の虚部

        rep = np.sqrt(complex(alpha,beta)) #比誘電率の平方根
        a = rep.real
        b = rep.imag
        proc = complex(-k0[i]*b,k0[i]*a) #伝搬定数

        kp = k0[i]*rep  #プラズマ中の波数
        a = kp.real
        b = kp.imag 
        ki = k0[i]      #真空中の波数
        cc = (ki+a)**2+b**2
        gamma = complex((ki**2-a**2-b**2)/cc,-2*ki*b/cc) #反射係数
        
        a = np.abs(gamma*(1-np.exp(-2.0*proc*d)))
        b = np.abs(1-gamma**2*np.exp(-2.0*proc*d))
        ref[i] = (a/b)  #反射率
        
        a = np.abs((1-gamma**2)*np.exp(-proc*d))
        b = np.abs(1-gamma**2*np.exp(-2.0*proc*d))
        trs[i] = (a/b)  #透過率
        
    ref = 20*np.log10(ref)
    trs = 20*np.log10(trs)
    return freq,ref,trs

def slab_plot():
    c  = 2.998e8      #光速
    freq = np.linspace(1e6,80e9,1000) #周波数
    w  = 2*np.pi*freq #角周波数
    nu = 20e9   #衝突周波数
    wp = 2*np.pi*28.7e9 #プラズマ周波数
    k0 = 2*np.pi*freq/c #波数
    d = 75e-6*200   #スラブ厚さ
    
    #ref   = np.zeros([len(freq)])   #反射強度
    #trs   = np.zeros([len(freq)])   #透過強度
    #nu = 0.0

    wc = 0
    freq,ref,trs = slab_th(wp,wc,nu,d)
    wc = 0.5
    freq,ref2,trs2 = slab_th(wp,wc,nu,d)
    plt.plot(freq,ref)
    plt.plot(freq,ref2)
    plt.vlines(28.7e9,-45,0,linestyle="dashed",color="red")
    plt.vlines(wc*28.7e9,-45,0,linestyle="dashed",color="blue")
    plt.ylim([-45,0])
    plt.xlim([0,80e9])
    plt.grid()
    plt.show()
    plt.close()


    plt.plot(freq,trs)
    plt.plot(freq,trs2)
    plt.vlines(28.7e9,-80,9,linestyle="dashed",color="red")
    plt.vlines(wc*28.7e9,-80,0,linestyle="dashed",color="blue")
    plt.ylim([-80,0])
    plt.xlim([0,80e9])
    plt.grid()
    plt.show()
    plt.close()

slab_plot()