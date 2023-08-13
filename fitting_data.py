import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize

number=np.array([40,60,80,100,120,140,160,180,200])
virial=np.array([0.706796,0.788093,0.869149,0.831645,0.975225,0.764741,0.900835,0.84093,0.873165])
error=np.array([0.000194908,0.000163204,0.000229463,0.000202798,0.000268487,0.000185878,0.000185878,0.000190708,0.000155824])

#define function for fit
def func(n,a,b):
    return a + b/np.sqrt(n)

sqrt_number_inv=1/np.sqrt(number)

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit,args=(sqrt_number_inv,virial,error), full_output=1)

pfinal = out[0]
covar = out[1]
print (pfinal)
print (covar)

index = pfinal[1]
amp = 10.0**pfinal[0]

indexErr = np.sqrt( covar[1][1] )
ampErr = np.sqrt( covar[0][0] ) * amp

##########
# Plotting data
##########

plt.clf()
plt.subplot(2, 1, 1)
plt.plot(number, virial)     # Fit
plt.errorbar(number, virial, yerr=error, fmt='k.')  # Data
plt.text(5, 6.5, 'Ampli = %5.2f +/- %5.2f' % (amp, ampErr))
plt.text(5, 5.5, 'Index = %5.2f +/- %5.2f' % (index, indexErr))
plt.title('Best Fit Power Law')
plt.xlabel('X')
plt.ylabel('Y')
#plt.xlim(1, 11)


# plt.subplot(2, 1, 2)
# plt.loglog(xdata, powerlaw(xdata, amp, index))
# plt.errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
# plt.xlabel('X (log scale)')
# plt.ylabel('Y (log scale)')
# plt.xlim(1.0, 11)

plt.show()