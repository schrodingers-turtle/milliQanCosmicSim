import random
import numpy as np
import scipy as sp
from scipy import integrate
import pylab as pl

#initialize cosmic variables
I_0   	= 88.0
n     	= 3
E_0   	= 4.29
E_c 	= 0.5
E_cmil	= 58
E_max 	= 4000.
epinv 	=  1/854.
Rod   	= 174.
N	= (n - 1)*(E_0 + E_c)**(n-1)

#initialize milliQan variables
depth 	= 60.	
srfrad	= 1000.		#surface radius for rmax
chmbrad = 2.		#chamber radius

#initialize simulation variables
rmin 	= 0
rmax 	= 1000
rstep	= 1
npart	= 1000

Imax	= 0.0085	#determined from Imax plots

def Dfun(theta):	#from ref[1]
	return (sp.sqrt(Rod**2*sp.cos(theta)**2 + 2*Rod + 1) - Rod*sp.cos(theta))**(-(n-1))

def Ifun(E,theta):	#also from ref[1]
	return I_0*N*(E_0 + E)**(-n)*(1 + E*epinv)**(-1)*Dfun(theta)

def test(I,E_min):	#for testing the sampled muons by energy distribution at theta = 0
	n_bins = 40
	Num = 100
	incr = sp.log10(E_max/E_min)/Num
	analytx = []
	analyty = []
	for i in range(0,Num+1):
		xx	= E_min*10**(i*incr)
		analytx.append(xx)

	analytx = np.asarray(analytx)

	
#	Norm = integrate.quad(I,E_min,E_max)[0]
	
#	analyty = np.asarray(analyty)/Norm

	j = 0

	samplx = []
	samply = []
	failx = []
	faily = []
	for i in range(1,npart):
		x = sp.log10(E_max/E_c)*random.random()
		y = Imax*random.random()
		xp = E_c*10**x
		yp = y
		if xp >= E_min and y <= I(xp):
			j += 1
			samplx.append(xp)
			samply.append(yp)
	
	pl.figure(figsize=(7,5))
	pl.plot(analytx, I(analytx), 'r', label='Analytical flux')
	pl.hist(samplx, n_bins, normed=1, label = 'Normed Muon count')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('E/GeV')
	ax.set_ylabel('I/s$^{-1}$(GeV/c)$^{-1}$')
	ax.legend()
#	pl.show()
	pl.savefig(str(rmin)+'to'+str(rmax)+'.jpg')
	pl.close()

def genmuons(I, E_min):
	Elist = []
	j = 0
	for i in range(0,npart*1000):
		x = sp.log10(E_max/E_c)*random.random()
		y = Imax*random.random()
		xp = E_c*10**x
		if y <= I(xp):
			j += 1
			if xp >= E_min:
				Elist.append([xp-E_min])
#	weight = integrate.quad(I,E_min,E_max)[0]/integrate.quad(I,E_c,E_max)[0]
	return [Elist,j]

def annul(r_a, r_b):
	r_mid 		= (r_a+r_b)/2.
	E_min		= sp.sqrt(r_mid**2+depth**2)-2				#min energy is 1 GeV per meter of rock
	theta 		= sp.arctan(r_mid/depth)
	
	area		= (r_b**2 - r_a**2)*np.pi

	thetap_b	= sp.arctan((r_mid - chmbrad)/depth)
	thetap_a	= theta + sp.arcsin(chmbrad/sp.sqrt(depth**2 + r_mid**2))
	phip		= sp.pi
	
	if r_mid > chmbrad:
		phip	= sp.arcsin(chmbrad/r_mid)

	sldangl 	= 2*phip*(sp.cos(thetap_b) - sp.cos(thetap_a))
	if thetap_b < 0:
		sldangl = 2*phip*(2 - sp.cos(thetap_b) - sp.cos(thetap_a))

#	print r_mid, E_min, theta, area, thetap_b, thetap_a, phip, sldangl

	def I(E):								#integrated flux for annulus r_ab
		return area*sldangl*Ifun(E,theta)
	
	rate = integrate.quad(I,E_min,E_max)[0]
	surfacerate = integrate.quad(I,E_c,E_max)[0]

#	test(I,E_min)								#produces a flux vs log momentum plot
	
	raw = genmuons(I,E_min)							#generates muons of varied momenta from normalized flux
	mlist = raw[0]
	allacc = raw[1]

	for i in range(0,len(mlist)):			
		phi = random.uniform(0,2*sp.pi)		
		mlist[i].append(random.uniform(r_a,r_b))			#gives muons random radial location within annulus
		mlist[i].append(phi)						#gives muons random azimuthal location
		mlist[i].append(random.uniform(thetap_a,thetap_b))		#gives muons random zenith trajectory within allowed range
		mlist[i].append((random.uniform(-phip,phip)+phi)%(2*sp.pi))	#gives muons random azimuthal trajectory within allowed range
#		mlist[i].append(I(E_cmil))					#calculates the maximum I at annulus to determine height of VN box
#		mlist[i].append(surfacerate)					#append surface rate to test for consistency
		mlist[i].append(rate)						#append rate to plot vs radius/azimuthal angle
#		mlist[i].append(weight)						#append weight per muon (obsolete)

	if 0==0:								#section for plotting histograms of generated muon energies
		grmlist = np.asarray(mlist)						
		Num = 100
		incr = sp.log10(E_max/E_cmil)/Num
		analytx = []
		for i in range(0,Num+1):
			xx	= E_cmil*10**(i*incr)
			analytx.append(xx)

		analytx = np.asarray(analytx)

		E 	= grmlist[:,0] + E_min
		w	= np.full(np.shape(E),rate/len(E))

		pl.figure(figsize=(7,5))
		pl.plot(analytx, I(analytx), 'r', label='Analytical flux')
		pl.hist(E, bins=100, weights=w, label = 'Normed Muon count')
		ax = pl.gca()
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_xlabel('E/GeV')
		ax.set_ylabel('I/s$^{-1}$(GeV)$^{-1}$')
		ax.legend()
	#	pl.show()
		pl.savefig(str(r_a)+'to'+str(r_b)+'.jpg')
		pl.close()

	return [mlist,allacc]							#returns final list of muons with radial and azimuthal location
										#and zenith and azimuthal trajectories
def main():

	def Dint(theta):							#Sanity check section
		return np.sin(theta)*Dfun(theta)
	
	angint = integrate.quad(Dint,0,np.pi/2)[0]

	def Iz(E):
		return Ifun(E,0)

	fullint = angint*integrate.quad(Iz,E_c,E_max)[0]*2*np.pi*0.0001*60	#rate per cm squared per minute

	print str(fullint) + " is the rate of muons on the surface per cm squared per minute."								

	Nsteps = (rmax - rmin)/rstep						#determine number of steps for looping over different annuli

	mlist = []
	allacc = 0

	while len(mlist) < npart:						#loop over annuli until desired number of muons is reached
		for i in range(1,Nsteps+1):					#loop over all annuli each iteration to prevent bias
			res = annul(rmin+(i-1)*rstep,rmin + i*rstep)
			mlist += res[0]
			allacc += res[1]

#	print(mlist)
	print str(len(mlist)/float(allacc)) + " is the fraction of muons that made it from the surface to milliQan."

#	np.save('muondat', np.asarray(mlist))



if __name__ == '__main__':
	main()

#references:
#[1] arXiv:1606.06907v3

