from pylab import *

def wbar(int,pl,atom):
    ''' determines the ionization rate for a pulse for a variety of noble gases. Uses the formulas in Zenghu Chang's book.

    input:
        int = laser intensity in 10**14 W cm**-2
        pl = pulse length in fs
        atom = "Xe", "Kr", etc.
    Output:
        cycle averaged ionization rate. '''

    F = sqrt(int/355)
    atom_dict = {"Xe": 0, "Kr": 1, "Ar": 2,"Ne":3}
    ind = atom_dict[atom]
    Ip = array([12.129, 13.999, 15.759, 21.564]);
    F0 = array([.84187, 1.04375, 1.24665, 1.99547])
    ns = array([1.0591, .9858, .92915, .7943])
    ls = array([0.0596, -0.01417, -0.07085, -0.2057])
    l = array([1, 1, 1, 1])
    m = array([0, 0, 0, 0])
    cnl = array([3.882, 4.025480, 4.11564, 4.2436])
    glm = array([3, 3, 3, 3])
    wb = sqrt(2/pi)*sqrt(3*F/2/F0)*cnl*glm*Ip*exp(-2*F0/3/F)*41.341*(2*F0/F)**(2*ns-m-1)
    return wb[ind]

def ionfrac(int, pl, atom):
    ''' determines the ionization fraction for a pulse for a variety of noble gases. Uses the formulas in Zenghu's book.

    input:
        int = laser intensity in 10**14 W cm**-2
        pl = pulse length in fs
        atom = "Xe", "Kr", etc.
    Output:
        cycle averaged ionization fraction. '''
    t = arange(-3*pl,3*pl,pl/10)
    lpulse = int*exp(-.5*(t/pl/2.355)**2)
    w_temp = zeros(size(t))
    j = 0
    for i in arange(0,size(t)):
        w_temp[i] = wbar(lpulse[i],pl,atom)
    w = sum(w_temp)*pl/10
    return 1-exp(-w)

def ioncrit(atom,q):
    ''' determines the critical ionization fraciton at 1070nm.
    input:
        atom = "Xe", "Kr", etc.
        q = harmonic order
    Output:
        critical ionization fraction at 1070 nm '''
    wr_dict = {"Xe": 10**15, "Kr": 2.195*10**15, "Ar": 2.83*10**15 , "Ne": 4.04*10**15}
    wr = wr_dict[atom]
    na_dict = {"Xe": 7.02*10**-4, "Kr": 4.27*10**-4, "Ar": 2.8*10**-4 , "Ne": .67*10**-4}
    na = na_dict[atom]
    lam = 1070e-9
    natm = 2.7*10**25
    re = 2.82*10**-15
    el = 1.6*10**-19
    me = 9.1*10**-31
    eps = 8.85*10**-12
    c = 3*10**8
    wp = sqrt(el**2*natm/eps/me/2)
    w = c*2*pi/(lam)
    dn = (na)-wp**2/(wr**2 - q**2*w**2)
    return (1+natm*re/2/pi/dn*lam**2)**-1
