import numpy as np
from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

plt.close('all')

w0 = 2*np.pi
ti = np.linspace(0,.25,200)
t = np.linspace(0,1,500)
tr = np.zeros(len(ti))

## Equations
def xfun(x, t0):
    return .5*(np.cos(w0*x)-np.cos(w0*t0)+w0*(x-t0)*np.sin(w0*t0))

def xfunanal(t):
    return (1/w0)*(np.pi/2-3*np.arcsin(2*w0*t/np.pi-1))

def KE(ti, tr):
    return 2*(np.sin(w0*tr)-np.sin(w0*ti))**2;

def laser(t):
    return np.cos(w0*t);

## Calculate return times
tr[0] =  fsolve(xfun, 1, args=(ti[0]))
for i in range(len(ti)-1):
    sol =  root(xfun, .9, args=(ti[i+1]))
    tr[i+1] = sol.x


## Plots
f, ax = plt.subplots()
ax.plot(ti,tr, 'r-', label = 'Numeric')
ax.plot(ti, xfunanal(ti), 'b-', label = 'Analytic')
ax.set_xlabel('Normalized Emission Time [arb.]')
ax.set_ylabel('Normalized Retrun Time [arb.]')
ax.legend(loc = 'lower left')
ax.set_ylim(0,1)
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
plt.show()

f, ax = plt.subplots()
ax.plot(ti,KE(ti,tr))
ax.set_xlabel('Emission Time')
ax.set_ylabel('Return Kinteic Energy [KE/U$_P$]')
ax.set_title('Semi-Classical Equation')
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
plt.show()

f, ax = plt.subplots()
ax.plot(KE(ti,tr),ti)
ax.set_xlabel('Emission Time')
ax.set_ylabel('Return Kinteic Energy [KE/U$_P$]')
ax.set_title('Semi-Classical Equation')
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
plt.show()

f, ax = plt.subplots()
ax.plot(t,laser(t),'k-', label = 'Laser')
ax.plot(np.linspace(.03,.77,100), xfun(np.linspace(.03,.77,100), .03), label = 'T$_0$ = .03')
ax.plot(np.linspace(.15,.45,100), xfun(np.linspace(.15,.45,100), .15), label = 'T$_0$ = .15')
ax.plot(np.linspace(.3,1,100), xfun(np.linspace(.3,1,100), .3), label = 'T$_0$ = .3')

ax.set_ylabel('Trajectory');
ax.set_xlabel('Time');
ax.set_title('Semi-Classical Equation');
plt.legend(loc = 'upper left')
plt.show()
