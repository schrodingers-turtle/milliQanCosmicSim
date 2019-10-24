import random
import scipy as sp
import pylab as pl

#initialize cosmic variables
I_0   	= 88.0
n     	= 3
E_0   	= 4.29
E_max 	= 4000.
epinv 	=  1/854.
Rod   	= 174.

#initialize milliQan variables
depth 	= 60.	
srfrad	= 1000.		#surface radius for rmax
chmbrad = 2.		#chamber radius

#initialize simulation variables
rmin 	= 0
rmax 	= 100
rstep	= 1
npart	= 1000

def Dfun(theta):	#from ref[1]
	return sp.sqrt(Rod**2*sp.cos(theta)**2 + 2*Rod + 1) - Rod*sp.cos(theta)

def Ifun(E,theta):	#also from ref[1]
	return I_0*(E_0 + E)**(-n)*(1 + E*epinv)**(-1)*Dfun(theta)**(-(n-1))

def test(I,E_min):	#for testing the sampled muons by energy distribution at theta = 0
	N = 100
	incr = sp.log10(E_max/E_min)/N
	analytx = []
	analyty = []
	for i in range(0,N+1):
		xx	= E_min*10**(i*incr)
		analytx.append(xx)
		analyty.append(I(xx))

	Imax = max(analyty)

	j		= 0

	samplx = []
	samply = []
	failx = []
	faily = []
	for i in range(1,npart*10):
		x = sp.log10(E_max/E_min)*random.random()
		y = Imax*random.random()
		xp = E_min*10**x
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

def genmuons(I, E_min):
	Imax = I(E_min)
	Elist = []
	for i in range(0,npart*10):
		x = sp.log10(E_max/E_min)*random.random()
		y = Imax*random.random()
		xp = E_min*10**x
		if y <= I(xp):
			Elist.append([xp])
	return Elist

def norm(I, E_min):
	return sp.quad(I,E_min,E_max)

def annul(r_a, r_b):
	r_mid 		= (r_a+r_b)/2.
	E_min		= sp.sqrt(r_mid**2+depth**2)-2	#min energy is 1 GeV per meter of rock
	theta 		= sp.arctan(r_mid/depth)
	
	area		= (r_b**2 - r_a**2)*sp.pi

	thetap_a	= sp.arctan((r_mid - chmbrad)/depth)
	thetap_b	= theta + sp.arcsin(chmbrad/sp.sqrt(depth**2 + r_mid**2))
	phip		= sp.pi
	if r_mid > chmbrad:
		phip	= sp.arcsin(chmbrad/r_mid)

	sldangl 	= 2*phip*(sp.cos(thetap_a) - sp.cos(thetap_b))

	def I(E):								#integrated flux for annulus r_ab
		return area*sldangl*Ifun(E,theta)
	
	def nI(E):								#normalized integrated flux for annulus r_ab
		return I(E)/norm(I,E_min)

#	test(nI,E_min)								#produces a flux vs log momentum plot
	
	mlist = genmuons(nI,E_min)						#generates muons of varied momenta from normalized flux
	
	for i in range(0,len(mlist)):			
		phi = random.uniform(0,2*sp.pi)		
		mlist[i].append(random.uniform(r_a,r_b))			#gives muons random radial location within annulus
		mlist[i].append(phi)						#gives muons random azimuthal location
		mlist[i].append(random.uniform(thetap_a,thetap_b))		#gives muons random zenith trajectory within allowed range
		mlist[i].append((random.uniform(-phip,phip)+phi)%(2*sp.pi))	#gives muons random azimuthal trajectory within allowed range
										#while adding azimuthal location and modding by 2pi

	return mlist								#returns final list of muons with radial and azimuthal location
										#and zenith and azimuthal trajectories
def main():
	Nsteps = (rmax - rmin)/rstep						#determine number of steps for looping over different annuli

	mlist = []

	while len(mlist) < npart/10:						#loop over annuli until desired number of muons is reached
		for i in range(1,Nsteps+1):					#loop over all annuli each iteration to prevent bias
			mlist += annul(rmin+(i-1)*rstep,rmin + i*rstep)
	print(mlist)

if __name__ == '__main__':
	main()

#references:
#[1] arXiv:1606.06907v3

