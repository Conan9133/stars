import pandas as pd
import numpy as np

dr21data=pd.read_csv("dr21all.csv")
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

#mean and standard deviation for vx and vy
mvx=-45.367184
svx=443.3824

mvy=-75.19972
svy=753.2898

#for selecting high velocity stars like within 2 sigma
flag=np.zeros(len(ra))

for i in range(len(ra)):
	if vx[i]> mvx-2*svx and vx[i]< mvx+2*svx :
		if vy[i]>mvy-2*svy and vy[i]< mvy+2*svy :
			flag[i]=1
			
 
dr21data.insert(loc=8,column="flag1",value=flag)
val={'vx':0,'vy':0}
dr21data=dr21data.fillna(value=val)

exportcsv=dr21data.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_dr21/dr21alldatanew.csv',index=None,header=True)

#find avg of vx and vy removing high velocity stars
s1=0
s2=0
c=0
for i in range(len(ra)):
	if flag[i]==1:
		s1=s1+vx[i]
		s2=s2+vy[i]
		c=c+1
avgvx=s1/c
avgvy=s2/c

print avgvx,avgvy,c,len(ra)-c

#counting matched stars with high velocity
vx1=dr21data['vx'];vx=vx1.values
count=0
for i in range(len(ra)):
	if vx[i]!=0:
		if flag[i]==0:
			count=count+1

print count

#counting stars with and without gaia matches
flgs=np.zeros(len(ra))
for i in range(len(ra)):
	if vx[i]!=0:
		flgs[i]=1
print flgs

with open('n.txt', 'w') as f:
    for item in flgs:
        f.write("%s\n" % item)

