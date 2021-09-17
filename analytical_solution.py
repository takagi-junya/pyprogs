from scipy.special import *
import matplotlib.pyplot as plt
import numpy as np
import phys_cons as pc

def cylinder_TE(k,a):
    m = 50
    dx = k[1]-k[0]
    k = k/a
    sigma = np.zeros([len(k)])
    
    u1 = np.zeros([2*m,len(k)])
    u2 = np.zeros([2*m,len(k)],dtype=np.complex128)
    u3 = np.zeros([2*m-1,len(k)-1])
    u4 = np.zeros([2*m-1,len(k)-1],dtype=np.complex128)
    
    for n in range(-m,m):
        for i in range(0,len(k)):
            u1[n,i] = jn(n,a*k[i])
            u2[n,i] = hankel2(n,a*k[i])
    
    for n in range(-m,m):
        u3[n,0] = (u1[n,1]-u1[n,0])/dx
        u4[n,0] = (u2[n,1]-u2[n,0])/dx
        for i in range(1,len(k)-1):
            u3[n,i] = (u1[n,i+1]-u1[n,i-1])/2*dx
            u4[n,i] = (u2[n,i+1]-u2[n,i-1])/2*dx
    
    for i in range(0,len(k)-1):
        sig = 0.0+0.0j
        for n in range(-m+1,m-1):
            c = (-1)**n
            u = u3[n,i]*u4[n,i].conjugate()
            l = u4[n,i]*u4[n,i].conjugate()
            sig = sig + c*u/l
        sigma[i]=((np.abs(sig)**2)*4/k[i])

    RCS = 10*np.log10(sigma/np.pi/a)
    return RCS

def cylinder_TM(k,a):
  m = 50
  k = k/a
  sigma = np.zeros([len(k)])
  for i in range(0,len(k)):
      sig = 0.0+0.0j
      for n in range(-m,m):
          u1 = (-1)**n
          u2 = jn(n,a*k[i])
          u3 = hankel2(n,a*k[i]).conjugate()
          l1 = np.abs(hankel2(n,a*k[i]))**2
          sig = sig+u1*u2*u3/l1
      sigma[i]=((np.abs(sig)**2)*4/k[i])
  RCS = 10*np.log10(sigma/np.pi/a)
  return RCS

def dhankel(n,x):
    return spherical_jn(n,x)-(0.0+1.0j)*spherical_yn(n,x)
    
def sphere_TE(ka):
    m=50
    x = ka
    sigma = np.zeros([len(ka)])
        
    for i in range(0,len(x)):
        sig = 0.0+0.0j
        b1  = 0.0+0.0j
        for n in range(1,m):
            a1 = (-1)**n
            a2 = (2*n+1)
            b1 = x[i]*dhankel(n,x[i])*(-n*dhankel(n,x[i])+x[i]*dhankel(n-1,x[i]))
            a  = a1*a2*b1.conjugate()
            b  = np.abs(b1)**2
            sig = sig + a/b
        sigma[i] = np.abs(sig)**2/(x[i]*x[i])

    RCS = 10*np.log10(sigma)
    RCS = sigma
    return RCS

def slab(wp,nu,d):
    c  = 2.998e8      #光速
    freq = np.linspace(1e6,100e9,1000) #周波数
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

def slab_mag_parallel(wp,wc,nu,d):
    c  = 2.998e8      #光速
    freq = np.linspace(1e6,100e9,1000) #周波数
    w  = 2*np.pi*freq #角周波数
    k0 = 2*np.pi*freq/c #波数
    
    ref   = np.zeros([len(freq)])   #反射強度
    trs   = np.zeros([len(freq)])   #透過強度

    print("wp:{:.3e}".format(wp))
    print("wc:{:.3e}".format(wc))
    print("nu:{:.3e}".format(nu))
    print("d:{:.3f}".format(d))
    
    for i in range(0,len(freq)):
        a = wp**2*(w[i]-wc)
        b = w[i]*((w[i]-wc)**2+nu**2)
        alpha = 1 - a/b  #比誘電率の実部
        beta  = -wp**2*nu/b   #比誘電率の虚部

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

def slab_mag_perpendicular(wp,wc,nu,d):
    c  = 2.998e8      #光速
    freq = np.linspace(1e6,100e9,1000) #周波数
    w  = 2*np.pi*freq #角周波数
    k0 = 2*np.pi*freq/c #波数
    
    ref   = np.zeros([len(freq)])   #反射強度
    trs   = np.zeros([len(freq)])   #透過強度

    print("wp:{:.3e}".format(wp))
    print("wc:{:.3e}".format(wc))
    print("nu:{:.3e}".format(nu))
    print("d:{:.3f}".format(d))
    
    for i in range(0,len(freq)):
        a = wp**2*((w[i]**2-wp**2)*(w[i]**2-wp**2-wc**2)+nu**2*w[i]**2)
        b = w[i]**2*(w[i]**2-wp**2-wc**2-nu**2)**2+nu**2*(2*w[i]**2-wp**2)**2
        alpha = 1 - a/b  #比誘電率の実部
        a = nu*wp**2*(wp**4+w[i]**2*(w[i]**2-2*wp**2+wc**2+nu**2))
        b = w[i]*(w[i]**2*(w[i]**2-wp**2-wc**2-nu**2)**2+nu**2*(2*w[i]**2-wp**2)**2)
        beta  = -a/b   #比誘電率の虚部

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
    freq = np.linspace(1e6,100e9,1000) #周波数
    w  = 2*np.pi*freq #角周波数k
    nu = 2*np.pi*3.18e9   #衝突周波数
    fp = 50e9 #プラズマ周波数
    wp = 2*np.pi*fp
    n = pc.eps0*pc.me*wp**2/(pc.qe**2)
    wc = 2*np.pi*47.75e9
    fc = wc/(2*np.pi)
    k0 = 2*np.pi*freq/c #波数
    d = (75e-6)*120   #スラブ厚さ
    
    #ref   = np.zeros([len(freq)])   #反射強度
    #trs   = np.zeros([len(freq)])   #透過強度
    #nu = 0.0

    freq,ref,trs = slab_mag_parallel(wp,wc,nu,d)
    freq,ref2,trs2 = slab_mag_perpendicular(wp,wc,nu,d)
    plt.plot(freq,ref,label="parallel")
    plt.plot(freq,ref2,label="perpendicular")
    plt.vlines(fp,-45,0,linestyle="dashed",color="red")
    plt.vlines(fc,-45,0,linestyle="dashed",color="blue")
    plt.ylim([-45,0])
    plt.xlim([0,100e9])
    plt.grid()
    plt.legend()
    plt.show()
    plt.close()


    plt.plot(freq,trs)
    plt.plot(freq,trs2)
    plt.vlines(fp,-80,9,linestyle="dashed",color="red")
    plt.vlines(fc,-80,0,linestyle="dashed",color="blue")
    plt.ylim([-80,0])
    plt.xlim([0,100e9])
    plt.grid()
    plt.show()
    plt.close()
    
#slab_plot()