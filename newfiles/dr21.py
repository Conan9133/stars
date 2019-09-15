import pandas as pd
import numpy as np

dr21data=pd.read_csv("dr21alldatanew.csv")

#counting stars with 0 vx or vy or pmra or pmdec
y=dr21data['pmra'];y1=y.values
d=0
for i in range(len(y1)):
	if  y1[i]==0:d=d+1
print d

x=dr21data['vx'];x1=x.values
c=0
for i in range(len(x1)):
	if  x1[i]==0:c=c+1
print c

#drawing a velocity from a normal distribution
dr21data.loc[dr21data.vx==0,"vx"]=np.random.normal(-45.367184,443.3824,c)
dr21data.loc[dr21data.vy==0,"vy"]=np.random.normal(-75.19972,753.2898,c)
dr21data.loc[dr21data.parallax==0,"parallax"]=1000/1500.0
dr21data.loc[dr21data.pmra==0,"pmra"]=np.random.normal(-0.5488128,13.341549,d)
dr21data.loc[dr21data.pmdec==0,"pmdec"]=np.random.normal(-2.285615,10.244071,d)

ra1=dr21data['ra'];ra=ra1.values
dec1=dr21data['dec'];dec=dec1.values
pmra1=dr21data['pmra'];pmra=pmra1.values
pmdec1=dr21data['pmdec'];pmdec=pmdec1.values
p1=dr21data['parallax'];p=p1.values
subc1=dr21data['subcluster'];subc=subc1.values
tar1=dr21data['target'];tar=tar1.values
rv1=dr21data['radial_velocity'];rv=rv1.values
vx1=dr21data['vx'];vx=vx1.values
vy1=dr21data['vy'];vy=vy1.values
flg1=dr21data['flag1'];flg=flg1.values
index1=dr21data.index
index=index1.to_numpy()

#counting stars with very high velocities
flag=np.zeros(len(ra1))
w=0
for i in range(len(ra1)):
	if vx[i]> mvx-2*svx and vx[i]< mvx+2*svx :
		if vy[i]>mvy-2*svy and vy[i]< mvy+2*svy :
			flag[i]=1
			w=w+1
 
dr21data.insert(loc=11,column="flag2",value=flag)
dr21data.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_dr21/dr21alldata1.csv',index=None,header=True)

