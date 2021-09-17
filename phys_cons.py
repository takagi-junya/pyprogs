import numpy as np
from scipy import constants as cons
c = cons.c
eps0 = cons.epsilon_0
mu0 = cons.mu_0
z0 = np.sqrt(mu0/eps0)
qe = cons.e
me = cons.m_e

def show():
    print("light speed:{:.3e}".format(c))
    print("eps0:{:.3e}".format(eps0))
    print("Z0:{:.3f}".format(z0))
    print("qe:{:.3e}".format(qe))
    print("me:{:.3e}".format(me))
