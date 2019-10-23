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

	def I(E):
		return area*sldangl*Ifun(E,theta)
	
#	test(I,E_min)		#produces a rate vs log momentum plot
	
	mlist = genmuons(I,E_min)
	
	for i in range(0,len(mlist)):
		phi = random.uniform(0,2*sp.pi)
		mlist[i].append(random.uniform(r_a,r_b))
		mlist[i].append(phi)
		mlist[i].append(random.uniform(thetap_a,thetap_b))
		mlist[i].append((random.uniform(-phip,phip)+phi)%(2*sp.pi))
	
	return mlist

def main():
	Nsteps = (rmax - rmin)/rstep

	mlist = []

	while len(mlist) < npart/10:
		for i in range(1,Nsteps+1):
			mlist += annul(rmin+(i-1)*rstep,rmin + i*rstep)
	print(mlist)

if __name__ == '__main__':
	main()

#references:
#[1] arXiv:1606.06907v3

