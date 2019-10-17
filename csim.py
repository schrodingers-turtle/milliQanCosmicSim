import random
import numpy as np
import pylab as pl

#initialize cosmic variables
I_0   	= 88.0
n     	= 3
E_0   	= 4.29
E_min	= 10. #cutoff value of Energy
E_max 	= 4000.
epinv 	=  1/854.
Rod   	= 174.
Norm	= (n-1)*(E_0 + E_min)**(n-1)

#initialize milliQan variables
depth 	= 60	
srfrad	= 1000	#surface radius for rmax
chmbrad = 2	#chamber radius

#initialize simulation variables
rmin 	= 0
rmax 	= 5
rstep	= 1
npart	= 10000

def Dfun(theta):
	return np.sqrt(Rod**2*np.cos(theta)**2 + 2*Rod + 1) - Rod*np.cos(theta)

def Ifun(E,theta):
	return I_0*Norm*(E_0 + E)**(-n)*(1 + E*epinv)**(-1)*Dfun(theta)**(-(n-1))

def test(I):	#for testing the sampled muons by energy distribution at theta = 0
	N = 101
	incr = (np.log10(4000)-1)/N
	analytx = []
	analyty = []
	for i in range(0,N):
		xx	= 10**(i*incr + 1)
		analytx.append(xx)
		analyty.append(I(xx))

	Imax = max(analyty)

	j		= 0
	Npts		= 1001

	samplx = []
	samply = []
	failx = []
	faily = []
	for i in range(1,Npts):
		x = (np.log10(4000)-1)*random.random()
		y = Imax*random.random()
		xp = 10**(x + 1)
		yp = y
		if y <= I(xp):
			j += 1
			samplx.append(xp)
			samply.append(yp)
		else:
			failx.append(xp)
			faily.append(yp)
	
	pl.figure(figsize=(7,5))
	pl.plot(analytx, analyty, 'b')
	pl.scatter(samplx, samply, c='g')
	pl.scatter(failx, faily, c='r')
	ax = pl.gca()
	ax.set_xscale('log')
	pl.show()

def genmuons(I):
	Imax = I(E_min)
	Elist = []
	for i in range(0,npart):
		x = (np.log10(4000)-1)*random.random()
		y = Imax*random.random()
		xp = 10**(x + 1)
		if y <= I(xp):
			Elist.append([xp])
	return Elist

def annul(r_a, r_b):
	r_mid 		= (r_a+r_b)/2.
	theta 		= np.arctan(r_mid/60.)
	phisub 		= 1.*(2*np.pi/360)
	theta_a		= np.arctan(r_a/60.)
	theta_b		= 2*theta - theta_a
	sldangl 	= phisub*(np.cos(theta_a) - np.cos(theta_b))
	area		= (r_b**2 - r_a**2)*phisub/2.

	thetap_a	= np.arctan((r_mid - 2.)/60.)
	thetap_b	= theta + np.arctan(2/np.sqrt(60**2 + r_mid**2))
	phip		= np.arctan(2/r_mid)

	def I(E):
		return area*sldangl*Ifun(E,theta)
	
	test(I)		#produces a rate vs log momentum plot
	
	mlist = genmuons(I)
	
	for i in range(0,len(mlist)):
		phi = random.uniform(0,2*np.pi)
		mlist[i].append(r_mid)
		mlist[i].append(phi)
		mlist[i].append(random.uniform(thetap_a,thetap_b))
		mlist[i].append(random.uniform(-phip,phip)+phi)
	
	return mlist

def main():
	Nsteps = (rmax - rmin)/rstep

	mlist = []

	while len(mlist) < npart:
		for i in range(1,Nsteps+1):
			mlist += annul(rmin+(i-1)*rstep,rmin + i*rstep)	
	print(mlist)

if __name__ == '__main__':
	main()
