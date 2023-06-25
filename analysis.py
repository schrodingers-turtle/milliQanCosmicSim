import numpy as np
import scipy as sp
import pylab as pl

n_bins = 100

mlist 	= np.load('muondat.npy')
E 	= mlist[:,0]
logE	= np.log10(E)
rate	= mlist[:,-1]
srfrate	= mlist[:,-2]
rad	= mlist[:,1]
sqrtrad	= np.sqrt(rad)
thetap	= mlist[:,3]
Imax	= mlist[:,-3]

print("there are " + str(mlist.shape[0]) + " muons")

totalrate = np.trapz(rate,x=rad)
totalsrfrt = np.trapz(srfrate,x=rad)

print("the total rate of muons hitting milliQan is " + str(totalrate) + " per second")

print("the total rate of muons at the surface is " + str(totalsrfrt) + " per second")

print("the ratio is " + str(totalrate/totalsrfrt))

pl.figure(figsize=(7,5))
pl.hist(E, n_bins)
ax = pl.gca()
ax.set_xlabel('E/GeV')
ax.set_ylabel('Number of muons')
#ax.legend()
#pl.show()
pl.savefig('summedhist.jpg')
pl.close()

pl.figure(figsize=(7,5))
pl.hist(E, n_bins)
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('E/GeV')
ax.set_ylabel('Number of muons')
#ax.legend()
#pl.show()
pl.savefig('summednumloghist.jpg')
pl.close()

pl.figure(figsize=(7,5))
n, b, patches = pl.hist(logE, n_bins)
bin_max = np.where(n == n.max())
print('maxbin', b[bin_max][0])
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('log of E/GeV')
ax.set_ylabel('Number of muons')
#ax.legend()
#pl.show()
pl.savefig('summedlogloghist.jpg')
pl.close()

pl.figure(figsize=(7,5))
pl.plot(rad,rate,',r')
ax = pl.gca()
ax.set_xlabel('r/m')
ax.set_ylabel('rate/s$^{-1}$')
#ax.legend()
#pl.show()
pl.savefig('ratevsrad.jpg')
pl.close()

pl.figure(figsize=(7,5))
pl.plot(thetap,rate,',g')
ax = pl.gca()
ax.set_xlabel('thetap/rad')
ax.set_ylabel('rate/s$^{-1}$')
#ax.legend()
#pl.show()
pl.savefig('ratevstheta.jpg')
pl.close()

pl.figure(figsize=(7,5))				#section for debugging and determining Imax for VN box height
pl.plot(rad,srfrate,',r')
ax = pl.gca()
ax.set_xlabel('r/m')
ax.set_ylabel('rate at surface/s$^{-1}$')
pl.savefig('srfratevsrad.jpg')
pl.close()

pl.figure(figsize=(7,5))
pl.plot(rad,Imax,',r')
ax = pl.gca()
ax.set_xlabel('r/m')
ax.set_ylabel('I$_{max}/s^{-1}$GeV$^{-1}$')
pl.savefig('IMaxvsrad.jpg')
pl.close()

#pl.figure(figsize=(7,5))
#pl.plot(rad,rate/srfrate,',r')
#ax = pl.gca()
#ax.set_xlabel('r/m')
#ax.set_ylabel('rate:srfrate ratio')
#pl.savefig('rateratiovsrad.jpg')
#pl.close()
