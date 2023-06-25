import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl

#initialize milliQan variables
depth 	= 60.	
srfrad	= 1000.		#surface radius for rmax
chmbrad = 2.		#chamber radius

#initialize muon mass
m_mu 	= 0.105658

mlist 	= np.load('muondat.npy')
E 	= mlist[:,0]
rad	= mlist[:,1]
phi	= mlist[:,2]
thetap	= mlist[:,3]
phip	= mlist[:,4]

E_min	= sp.sqrt(rad**2+depth**2)-2
srfE	= E + E_min

skinvar	= np.zeros((np.shape(mlist)[0],7))	#kinetic variables at surface  (x,y,z,E,px,py,pz)	
mkinvar	= np.zeros((np.shape(mlist)[0],7))	#kinetic variables at milliQan (x,y,z,E,px,py,pz)

#start in frame where phi = 0 then rotate to final desired phi value
phip0	= phip - phi

#set unrotated surface position variables
skinvar[:,0] = rad
skinvar[:,1] = 0
skinvar[:,2] = depth
skinvar[:,3] = srfE

#set slopes for solver
zxslope	= 1/np.tan(thetap)
yxslope = -np.tan(phip0)/np.sin(thetap)

#solve for points on sphere of milliQan
for i in range(0,np.shape(mlist)[0]):
	slopesq = zxslope[i]**2+yxslope[i]**2
	mx 	= np.roots([1+slopesq,2*depth*zxslope[i]-2*rad[i]*slopesq,slopesq*rad[i]**2-2*depth*zxslope[i]*rad[i]+depth**2-chmbrad**2])
	if isinstance(mx[0],complex):			#filter out particles that miss
		mkinvar[i,0] = np.NaN
	elif zxslope[i]*(mx[0] - rad[i]) >= -depth and zxslope[i]*(mx[0] - rad[i]) > zxslope[i]*(mx[1] - rad[i]):	#ensure all strikes are on top hemisphere and entry point is logged
		mkinvar[i,0] = mx[0]
	elif zxslope[i]*(mx[1] - rad[i]) >= -depth:	
		mkinvar[i,0] = mx[1]
	else:
		mkinvar[i,0] = np.NaN			#get rid of all other points
			
mkinvar[:,1] = yxslope*(mkinvar[:,0] - rad)				#set milliQan position variables
mkinvar[:,2] = zxslope*(mkinvar[:,0] - rad) + depth

mkinvart     = mkinvar[:,0]*np.cos(phi) - mkinvar[:,1]*np.sin(phi)	#rotate milliQan positions to phi frame
mkinvar[:,1] = mkinvar[:,0]*np.sin(phi) + mkinvar[:,1]*np.cos(phi)
mkinvar[:,0] = mkinvart

mkinvar[:,3] = E							#set energy

skinvart     = skinvar[:,0]*np.cos(phi) - skinvar[:,1]*np.sin(phi)	#rotate surface positions to phi frame
skinvar[:,1] = skinvar[:,0]*np.sin(phi) + skinvar[:,1]*np.cos(phi)
skinvar[:,0] = skinvart

mkinvar[:,4] = mkinvar[:,0] - skinvar[:,0]				#set preliminary vectors based on displacement from surface to milliQan sphere
mkinvar[:,5] = mkinvar[:,1] - skinvar[:,1]
mkinvar[:,6] = mkinvar[:,2] - skinvar[:,2]

veclen = np.sqrt(mkinvar[:,4]**2 + mkinvar[:,5]**2 + mkinvar[:,6]**2)

mkinvar[:,4] = mkinvar[:,4]/veclen
mkinvar[:,5] = mkinvar[:,5]/veclen
mkinvar[:,6] = mkinvar[:,6]/veclen

smom = np.sqrt(srfE**2 + 2*srfE*m_mu)

skinvar[:,4] = mkinvar[:,4]*smom
skinvar[:,5] = mkinvar[:,5]*smom
skinvar[:,6] = mkinvar[:,6]*smom

mmom = np.sqrt(E**2 + 2*E*m_mu)

mkinvar[:,4] = mkinvar[:,4]*mmom
mkinvar[:,5] = mkinvar[:,5]*mmom
mkinvar[:,6] = mkinvar[:,6]*mmom

#print mkinvar[:,0]**2 + mkinvar[:,1]**2 + mkinvar[:,2]**2

print(mkinvar)
print(skinvar)

rate	= mlist[:,5]
uqind	= np.unique(rate,return_index=True)[1]


fig = pl.figure()							#3d plots for debugging
ax = fig.add_subplot(111, projection='3d')
ax.scatter(mkinvar[:,0], mkinvar[:,1], mkinvar[:,2], marker='.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#pl.show()
pl.savefig('3dscatter.jpg')
pl.close()

fig = pl.figure()							#3d plots for debugging
ax = fig.add_subplot(111, projection='3d')
ax.quiver(mkinvar[:,0][uqind], mkinvar[:,1][uqind], mkinvar[:,2][uqind], mkinvar[:,4][uqind], mkinvar[:,5][uqind], mkinvar[:,6][uqind], length=0.1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#pl.show()
pl.savefig('3dvec.jpg')
pl.close()

np.save('mQvecs', np.asarray(mkinvar))
np.save('srfvecs', np.asarray(skinvar))
