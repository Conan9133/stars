import pandas as pd
import numpy as np
import numpy.random as rnd
from numpy import pi, sqrt

dr21data=pd.read_csv("dr21Inew.csv")
ra1=dr21data['ra'];rastars=ra1.values
dec1=dr21data['dec'];decstars=dec1.values
vxavg=-13.47146225584369
vyavg=-28.75265571764528
vx1=dr21data['vx'];vx=vx1.values
vy1=dr21data['vy'];vy=vy1.values
vx=vx-vxavg
vy=vy-vyavg
rastars=rastars*pi/180.0
decstars=decstars*pi/180.0
nobsstars=len(rastars)
dcluster=1500.0
racentre=309.7536344*pi/180.0  
deccentre=42.4118879*pi/180.0
rcminor=0.25*0.25*60/206265*dcluster
obsstars_zfactor=4

cluster_ra_avg=309.752083*pi/180.0
cluster_dec_avg=42.345092*pi/180.0
cluster_zmax=0.0

dra     = rastars - racentre
cosdec  = np.cos(deccentre)
sindec  = np.sin(deccentre) 
x       = (np.cos(decstars) * np.sin(dra) /(cosdec * np.cos(decstars) * np.cos(dra) + sindec * np.sin(decstars)))
y       = ((cosdec * np.sin(decstars) - sindec * np.cos(decstars) * np.cos(dra))/(sindec * np.sin(decstars) +cosdec * np.cos(decstars) * np.cos(dra)))      
newx = x * dcluster
newy= y * dcluster
newz= rnd.uniform(-obsstars_zfactor * rcminor,obsstars_zfactor  * rcminor,nobsstars)

cosdec  = np.cos(cluster_dec_avg)
sindec  = np.sin(cluster_dec_avg)
dra     = racentre -cluster_ra_avg
dx      = (np.cos(deccentre) * np.sin(dra) /(cosdec * np.cos(deccentre) * np.cos(dra) +sindec * np.sin(deccentre)))
dy      = ((cosdec * np.sin(deccentre) -np.cos(deccentre) * sindec * np.cos(dra))/(sindec * np.sin(deccentre) +cosdec * np.cos(deccentre) * np.cos(dra)))
dx     *= dcluster
dy     *= dcluster
dz      = rnd.uniform(-1.0*cluster_zmax,cluster_zmax)


dr21data['newx']=newx+dx
dr21data['newy']=newy+dy
dr21data['newz']=newz+dz


exportcsv=dr21data.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_dr21/newdata/dr21Inew1.csv',index=None,header=True)







