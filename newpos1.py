import pandas as pd
import numpy as np
import numpy.random as rnd
from numpy import pi, sqrt

dr21data=pd.read_csv("dr21Xnew.csv")
ra1=dr21data['ra'];rastars=ra1.values
dec1=dr21data['dec'];decstars=dec1.values
rastars=rastars*pi/180.0
decstars=decstars*pi/180.0
nobsstars=len(rastars)
dcluster=1500.0
#racentre=309.7723832*pi/180.0  
#deccentre=42.3549951*pi/180.0
#rcminor=0.80*0.25*60/206265*dcluster
obsstars_zfactor=4

cluster_ra_avg=309.752083*pi/180.0
cluster_dec_avg=42.345092*pi/180.0
cluster_zmax=0.0
unclustered_zdistance=2.0

dra     = rastars - cluster_ra_avg
cosdec  = np.cos(cluster_dec_avg)
sindec  = np.sin(cluster_dec_avg)
x = (np.cos(decstars) * np.sin(dra) /(cosdec * np.cos(decstars) * np.cos(dra) + sindec * np.sin(decstars)))
y = ((cosdec * np.sin(decstars) -sindec * np.cos(decstars) * np.cos(dra))/(sindec * np.sin(decstars) +cosdec * np.cos(decstars) * np.cos(dra)))
newx = x * dcluster 
newy = y * dcluster 
newz = rnd.uniform(-1.0*unclustered_zdistance,unclustered_zdistance,nobsstars)

dr21data['newx']=newx
dr21data['newy']=newy
dr21data['newz']=newz
exportcsv=dr21data.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_dr21/newdata/dr21Xnew1.csv',index=None,header=True)

