# QuantumPaths.py
from pylab import *
from scipy.optimize import root, fsolve
from matplotlib.ticker import LinearLocator
#This script solves the three equations obtained from minimizing the
#quasiclassical action. See Mairesse_Science2003 supplementary material.
#Craig Benko

close('all')


E = np.linspace(0.01,.07,50)
om = 1.2/27.2
Ip = 12.13/27.2 # IP of Xe in atomic units
q = np.arange(13,31,2)
qplotme = 19
check = 0

ps = np.zeros((len(q),len(E)), dtype = 'complex128')
tis = np.zeros((len(q),len(E)), dtype = 'complex128')
trs = np.zeros((len(q),len(E)), dtype = 'complex128')

pl = np.zeros((len(q),len(E)), dtype = 'complex128')
til = np.zeros((len(q),len(E)), dtype = 'complex128')
trl = np.zeros((len(q),len(E)), dtype = 'complex128')

def saddle(X,E0,om,Ip,q):
    p  = X[0]
    ti = X[1]
    tr = X[2]

    Ati = E0*np.sin(om*ti)/om
    Atr = E0*np.sin(om*tr)/om

    Y1 = 0.5*(p-Ati)**2+Ip
    Y2 = p*(tr-ti) + E0/om**2 *  (np.cos(om*tr)-np.cos(om*ti))
    Y3 = 0.5*(p-Atr)**2 + Ip - q*om

    Y = [Y1,Y2,Y3]
    return Y

def real_saddle(X,E0,om,Ip,q):
    XX = [X[0] + 1j*X[1], X[2]+ 1j*X[3], X[4]+1j*X[5] ]
    real_fun = saddle(XX,E0,om,Ip,q)
    return [real(real_fun[0]),imag(real_fun[0]),real(real_fun[1]),imag(real_fun[1]), real(real_fun[2]),imag(real_fun[2])]

def cutoff(E0, om, Ip):
    # calculate classical cutoff energy
    Ec = 3.17 * E0**2 / (4 * om**2) + Ip
    qc = floor(Ec / om)
    return Ec, qc

def plateau(q, om, Ip):
    return ((q*om-Ip)*4*om**2/3.17)*350

# calculate classical cutoff energy
Ec, qc = cutoff(E, om, Ip)

# start looping
for k in range(len(E)):
    E0 = E[k]

    X0s = .95*array([1.5, 0,  17.2, 20.0 , 60.5, -1]) # short guess
    X0l = .9*array([0.1, 0, 0.87, 15.6, 132.2, 0.1]) # long guess
    Xs = np.zeros((len(q),6)) # Short Traj
    Xl = np.zeros((len(q),6)) # Long Traj

    Xs[0,:] = root(real_saddle, X0s, args=(E0,om,Ip,q[0]), method = 'lm', tol = 1e-4).x
    Xl[0,:] = root(real_saddle, X0l, args=(E0,om,Ip,q[0]), method = 'lm', tol = 1e-4).x

    for kk  in range(len(q)-1):
        Xsguess = Xs[kk,:]
        Xlguess = Xl[kk,:]
        Xs[kk+1,:]  = root(real_saddle, X0s, args=(E0,om,Ip,q[kk+1]), method = 'hybr', tol = 1e-4).x
        Xl[kk+1,:] = root(real_saddle, X0l, args=(E0,om,Ip,q[kk+1]), method = 'hybr', tol = 1e-4).x

    # unpack answers
    ps[:,k]  = Xs[:,0] + Xs[:,1]*1j
    tis[:,k] = Xs[:,2] + Xs[:,3]*1j
    trs[:,k] = Xs[:,4] + Xs[:,5]*1j

    pl[:,k]  = Xl[:,0] + Xl[:,1]*1j
    til[:,k] = Xl[:,2] + Xl[:,3]*1j
    trl[:,k] = Xl[:,4] + Xl[:,5]*1j

if check == 1:
    f, ax = plt.subplots()
    ax.plot(q,real(trs[:,10])*24.2,'k', label = 'short')
    ax.plot(q,real(trl[:,10])*24.2,'r', label = 'long')
    axvline(qc[10])
    plt.legend()
    ax.set_title(str(E[10]**2*350) + 'x10$^{14}$ W/cm$^2$')
    ax.set_xlabel('Harmonic Order');
    ax.set_ylabel('Emission Time [as]');
    plt.show()

    f, ax = plt.subplots()
    ax.plot(q,real(trs[:,20])*24.2,'k', label = 'short')
    ax.plot(q,real(trl[:,20])*24.2,'r', label = 'long')
    axvline(qc[20])
    plt.legend()
    ax.set_title(str(E[20]**2*350) + 'x10$^{14}$ W/cm$^2$')
    ax.set_xlabel('Harmonic Order');
    ax.set_ylabel('Emission Time [as]');
    plt.show()
else:
    pass

taus = trs-tis # time in continuum! Note that it is complex.
taul = trl-tis

Up = E**2/(4*om**2)
IWcm2 = 1e-4*1/2*1/377*5.142e11**2*E**2

# sort trajectories
q = qplotme
loc = where(q == qplotme)[0]


S1s = q*om*trs[loc,:]
S2s = -Ip*taus[loc,:]
S3s = -.5*ps[loc,:]**2*taus[loc,:]
S4s = ps[loc,:]*E/om**2*(np.cos(om*tis[loc,:])-np.cos(om*trs[loc,:]))
S5s = -.25*E**2/om**2*taus[loc,:]
S6s = E**2/(8*om**3)*(np.sin(2*om*trs[loc,:])-sin(2*om*tis[loc,:]))
Ss = S1s+S2s+S3s+S4s+S5s+S6s
phis = real(Ss)[0]
phis = unwrap(phis-phis[0],pi)
# long trajectories
S1l = q*om*trl[loc,:]
S2l = -Ip*taul[loc,:]
S3l = -.5*pl[loc,:]**2*taul[loc,:]
S4l = pl[loc,:]*E/om**2*(np.cos(om*til[loc,:])-np.cos(om*trl[loc,:]))
S5l = -.25*E**2/om**2*taul[loc,:]
S6l = E**2/(8*om**3)*(np.sin(2*om*trl[loc,:])-np.sin(2*om*til[loc,:]))
Sl = S1l+S2l+S3l+S4l+S5l+S6l
phil = real(Sl)[0]
phil = unwrap(phil-phil[0],pi)

f, ax = plt.subplots()

ax.plot(Up*4*om**2*35000,phis,'k', label = 'Short')
ax.plot(Up*4*om**2*35000,phil,'r', label = 'Long')
ax.plot(Up*4*om**2*35000,phis+phil,'b', label = 'Total')
axvline(plateau(qplotme,om, Ip)*100, color = 'm', label = 'Cutoff Intensity')

ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
ax.set_xlabel('Intensity [TW/cm$^2$]')
ax.set_ylabel(r'$\phi$ [rad]')
ax.legend(loc ='lower center')
ax.set_title(str(q) + '$^{th}$ Harmonic')
plt.show()


f, ax = plt.subplots()
ax.plot(Up[:-1]*4*om**2*35000,diff(phis)/diff(Up*4*om**2*35000),'k', label = 'Short')
ax.plot(Up[:-1]*4*om**2*35000,diff(phil)/diff(Up*4*om**2*35000),'r', label = 'Long')
axvline(plateau(qplotme,om, Ip)*100, color = 'm', label = 'Cutoff Intensity')

ax.set_ylim(-.5,.5)
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
ax.set_xlabel('Intensity [TW/cm$^2$]')
ax.set_ylabel(r'$\alpha$ [rad TW$^{-1}$ cm$^2$]')
ax.legend(loc ='upper right')
ax.set_title(str(q) + '$^{th}$ Harmonic')
plt.show()
