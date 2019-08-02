#Consider MYStIX as leading catalogue & all the angles are in degrees

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import csv
from astropy.io import ascii
import warnings
import math
from collections import Counter
a=0.000277778

#extracting cluster/region data from whole data of MYStIX
target=input("Name of cluster/region:")
#dat=Table.read('apjs485825_MPCM.fits',format='fits')
data=pd.read_csv('mystixfull.csv')
cluster=data[data['target']==target]
cluster1=pd.DataFrame(cluster)

print len(cluster1)

#opening MYStIX data & sorting according to declination
#mystixsort=cluster1.sort_values(['dec'],ascending=False)
#mystixsort=mystixsort.reset_index(drop=True)
values1={'MAG_J':0,'Class_Pos_Err':0}
cluster1=cluster1.fillna(value=values1)
items=['99','-99','2MASS_G05','Megeath2012','2MASS_MC04','VLT','MLLA_NTT','MLLA_FLWO','OBcat','DBB','AEE','OV','2011ApJS..194...10P','UAA','0cc','ABB','AAA','TWOMASS','HC2000','OO','ON','UAU','CBA','CAB','Spitzer Vela','...']
cluster1.loc[cluster1['Class_Pos_Err'].isin(items),'Class_Pos_Err']=0
cluster1.loc[cluster1['MAG_J'].isin(items),'MAG_J']=0
mystixsortra=cluster1['ra'];mystixra=mystixsortra.values
mystixsortdec=cluster1['dec'];mystixdec=mystixsortdec.values
mystixsortmag=cluster1['MAG_J'];mystixjmag=mystixsortmag.values
index1=cluster1.index
indexmys=index1.to_numpy()

#finding maximum error in mystix data

mystixsorterr=cluster1['Class_Pos_Err'];mystixerr=mystixsorterr.values
export_csv=cluster1.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/cluster.csv',index=None,header=True)

#finding the range(check later whether the range is reasonable or not)??????////////////////////
p=float(input("Enter excess box length:"))
minra=min(mystixra)-p
maxra=max(mystixra)+p
rra=maxra-minra
mindec=min(mystixdec)-p
maxdec=max(mystixdec)+p
rdec=maxdec-mindec
print rra,rdec

#accesing gaia data from archive
target_name=input("Enter the exact name:")
clusterdata=SkyCoord.from_name(target_name)
ra1=clusterdata.ra
dec1=clusterdata.dec
coord=SkyCoord(ra=ra1,dec=dec1,unit=(u.degree,u.degree),frame='icrs')
width = u.Quantity(rra, u.deg)
height = u.Quantity(rdec, u.deg)
r=Gaia.query_object_async(coordinate=coord, width=width, height=height)
ascii.write(r,'gaia.csv',exclude_names=['designation'],delimiter=',')

values={'pm':0.00005}
#opening gaia data & sorting data according to declination & adding pm column
gaia=pd.read_csv('gaia.csv',usecols=['source_id','ra','ra_error','dec','dec_error','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error','ra_dec_corr','phot_g_mean_mag','phot_rp_mean_mag','a_g_val','astrometric_excess_noise','bp_rp','radial_velocity','radial_velocity_error','phot_variable_flag','teff_val'])
w1=gaia['ra'];w=w1.values
ww=len(w)
gaia=gaia.dropna(subset=['ra','dec','pmra','pmdec'])
gaia=gaia.fillna(value=values)
gaiasort=gaia.sort_values(['dec'],ascending=False) #sorting
gaiasort=gaiasort.reset_index(drop=True) #reseting index
x=gaiasort.apply(lambda row:math.sqrt((row.pmra)**2+(row.pmdec)**2),axis=1) #pm
val=x.values #accesing pm values
gaiasort.insert(loc=11,column='pm',value=val) #inserting pm to dataframe
#gaiasort=gaiasort.loc[gaiasort['astrometric_excess_noise']<1]
#gaiasort=gaiasort.reset_index(drop=True)
gaiasortra=gaiasort['ra'];gaiara=gaiasortra.values 
gaiasortdec=gaiasort['dec'];gaiadec=gaiasortdec.values
raerr1=gaiasort['ra_error'];raerr=raerr1.values
decerr1=gaiasort['dec_error'];decerr=decerr1.values
radeccorr1=gaiasort['ra_dec_corr'];radeccorr=radeccorr1.values
pm1=gaiasort['pm']
pm=pm1.values
gaiamag1=gaiasort['phot_rp_mean_mag'];gaiarmag=gaiamag1.values
index2=gaiasort.index
indexgaia=index2.to_numpy()

print len(gaiara)

#determine position errors
#mme=1.035*a #maximum mystix error in arcsec
PosErr=[]
for i in range(len(mystixra)):
	err=float(mystixerr[i])*a
	raerrmax=max(raerr)
	decerrmax=max(decerr)
	maxm=max(raerrmax,decerrmax)
	PosErr.append(err+maxm*a*0.001)
print len(PosErr)
print max(PosErr)

#determine each first filter radius for each object
#Rden=30*a #arcsec
radius=[]
H=3
EpochDiff=12 #years
for i in range(len(mystixra)):
	pmmax=max(pm)
	w=H*PosErr[i] + (pmmax*a*EpochDiff*10e-4)
	radius.append(w)	
print len(radius)
print max(radius)

#finding neighbours in bigger circle (Candidates)
getindex=[]
for i in range(len(mystixra)):	#mystixdata
	getindex.append([])
	r=radius[i]
	d=0
	for j in range(len(gaiara)):	#gaiadata
		if mystixdec[i]-r<gaiadec[j] and gaiadec[j]<mystixdec[i]+r :
			rag=gaiara[j]	
			decg=gaiadec[j]
			ram=mystixra[i]
			decm=mystixdec[i]
			d=2*math.asin((math.sin((decm-decg)*0.5))**2+((math.sin((ram-rag)*0.5))**2)*(math.cos(decg))*(math.cos(decm)))
			if d < r:
				getindex[i].append(indexgaia[j])
			else:
				continue
		elif mystixdec[i]-r>gaiadec[j]:
			break

print len(getindex)

#finding local density for each mystix object
dens=[]
for i in range(len(getindex)): #len==number of mystix objects
	c=0
	for j in range(len(getindex[i])):
		c+=1
	density=c/(math.pi*radius[i]**2)
	dens.append(c)

#print dens
print min(dens),max(dens)
#print len(dens)

#members inside small circle (Good neighbours) 
pmeff=0.2*a	#threshold proper motion
goodindex=[]
for i in range(len(getindex)):	#mystix
	goodindex.append([])
	for j in range(len(getindex[i])):	#gaia
		ram1=mystixra[i]	
		decm1=mystixdec[i]		
		rag1=gaiara[getindex[i][j]]
		decg1=gaiadec[getindex[i][j]]
		ragerr=raerr[j]*0.001*a
		pmg=pm[j]*0.001*a	
		gaiaerr=ragerr+pmg*EpochDiff*0.2		
		#for mystix error
		#sysmys=pmeff*EpochDiff*0.2
		#myserr=math.sqrt(mme**2+sysmys**2)
		myserr1=float(mystixerr[i])*a*0.707
		toterr=math.sqrt(gaiaerr**2+myserr1**2) #total error
		#for distance
		dist=2*math.asin((math.sin((decm1-decg1)*0.5))**2+((math.sin((ram1-rag1)*0.5))**2)*(math.cos(decg1))*(math.cos(decm1)))
		distrat=dist/toterr
		if distrat<=3.43:
			goodindex[i].append(getindex[i][j])
#			print toterr,dist,distrat

print len(goodindex)

#Figure of merit(FOM)
score=[]
for i in range(len(goodindex)):	#mystix
	ram2=mystixra[i]	
	decm2=mystixdec[i]
	locdens=dens[i]
	score.append([])
	for j in range(len(goodindex[i])):	#gaia
		rag2=gaiara[goodindex[i][j]]
		decg2=gaiadec[goodindex[i][j]]		
		sxm=float(mystixerr[i])*a*0.707
		sym=float(mystixerr[i])*a*0.707
		sxg=raerr[j]*a*0.001	
		syg=decerr[j]*a*0.001	
		sxc2=sxm**2+sxg**2
		syc2=sym**2+syg**2
		sm=sxc2
		distance=2*math.asin((math.sin((decm2-decg2)*0.5))**2+((math.sin((ram2-rag2)*0.5))**2)*(math.cos(decg2))*(math.cos(decm2)))
		fr=distance/(sxc2*1.0*radeccorr[j])
		fom=(math.exp(-0.5*fr**2))/(2*math.pi*locdens)		 
		sc=math.asinh(fom)
		score[i].append(sc)
		#print "score[",i,"][",j,"]=",score[i][j]

#mystix objects having highest score
w2=0
best=[]
scoremax=[]
for i in range(len(score)):
	if len(score[i])!=0:	
		scoremax.append(max(score[i]))		
		ind_max=np.argmax(score[i])	
		best.append(goodindex[i][ind_max])			
		w2+=1
	else:
		scoremax.append(0)
		best.append(0)

#print w2
#print best #index of gaia members
#for i in range(len(best)):
#	print "best[",i,"]=",best[i]
#	print "scoremax[",i,"]=",scoremax[i]

print len(best)
print max(best),min(best)

#finding possible duplicates in MYStIX
dpm=[]
dpg=[]
best1=set(best)
best2=list(best1)
#best2.remove(0)
best1=best2
for k in range(len(best1)):
	dpm.append([])
	dpg.append([])
	for i,j in enumerate(best):
		if j==best1[k]:
			dpm[k].append(i)
			dpg[k].append(j)

#print dpm  #indices of MYStIX
print "---------------------------------------\n"
#print dpg  #indices of GAIA
print len(dpm),len(dpg)

#File for particular region//////////////////////////////////////////////////
#create dataframe & appending rows
df=pd.DataFrame(columns=['source_id','ra','ra_error','dec','dec_error','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error','phot_g_mean_mag','phot_rp_mean_mag','a_g_val','astrometric_excess_noise','bp_rp','radial_velocity','radial_velocity_error','phot_variable_flag','teff_val'])

df1=pd.DataFrame(columns=['designation','ra','dec','xray','ir','ob','me','color','lum','target','subcluster','MYSTIX_SFR','Class_Name','Class_RAdeg','Class_DEdeg','Class_Pos_Err','Class_Pos_Origin','XRAY_NAME','XRAY_LABEL','NIR_NAME','NIR_LABEL','MIR_NAME','MIR_LABEL','OB3_LABEL','XCAT_INDEX','ISED_INDEX','SPTY','ORIGIN_SpTy','MAG_OB','BAND_OB','ProbNoSrc_min','ProbKS_single','ProbKS_merge','ExposureTimeNominal','ExposureFraction','NumObservations','NumMerged','Theta_Lo','Theta','Theta_Hi','PsfFraction','AfterglowFraction','SrcCounts_t','NetCounts_t','NetCounts_s','NetCounts_h','NetCounts_Lo_t','NetCounts_Hi_t','NetCounts_Lo_s','NetCounts_Hi_s','NetCounts_Lo_h','NetCounts_Hi_h','MedianEnergy_t','log_PhotonFlux_t','log_PhotonFlux_s','log_PhotonFlux_h','LX_H','LX_HC','SLX_HC_STAT','SLX_HC_SYST','LX_T','LX_TC','SLX_TC_STAT','SLX_TC_SYST','LOGNH_OUT','SLOGNH_OUT_STAT_OUT','SLOGNH_OUT_SYST_OUT','XN_PROB_CP','XM_PROB_CP','MAG_J','ERROR_J','MAG_H','ERROR_H','MAG_K','ERROR_K','MAG_3p6um','ERROR_3p6um','MAG_4p5um','ERROR_4p5um','MAG_5p8um','ERROR_5p8um','MAG_8p0um','ERROR_8p0um','J_FLAG','H_FLAG','K_FLAG','CC_FLG','PH_QUAL','SQF_J','SQF_H','SQF_K','SQF_3P6UM','SQF_4P5UM','SQF_4P8UM','SQF_8P0UM','AP_LS_FLG','ORIGIN_J','ORIGIN_H','ORIGIN_K','ORIGIN_3p6um','ORIGIN_4p5um','ORIGIN_5p8um','ORIGIN_8p0um','SED_FLG','SED_AV','SED_STAGE','H1_prior','H2_prior','H3_prior','H4_prior','H1_posterior','H2_posterior','H3_posterior','H4_posterior','H2_dominant_factor','xray_class_code'])
#accessing members
count=0
for i in range(len(dpm)):
	if len(dpm[i])==0:
		continue
	elif len(dpm[i])==1:
		x=dpg[i][0]
		y=dpm[i][0]
		df=df.append(gaiasort.loc[x,:])
		df1=df1.append(cluster1.iloc[y,:])
		#print "gaia[",i,"][0]=",dpg[i][0]
		#print "mystix[",i,"][0]=",dpm[i][0]  
	elif len(dpm[i])!=1:
		sml=[]
		x1=[]
		for j in range(len(dpm[i])):
			x=dpm[i][j]		#accesing mystix index
			sm1=scoremax[x]		#finding score of mystix data
			sml.append(sm1)		#adding it to array
			x1.append(x)		#storing mystix indices
		smax=max(sml)			#finding maximum score among mystix objects
		if smax==0:
			magtot=[]
			for k in range(len(x1)):
				if mystixjmag[x1[k]]!=0 or gaiarmag[dpg[i][k]]!=0:
					mag1=float(mystixjmag[x1[k]])-gaiarmag[dpg[i][k]]
					magtot.append(mag1)
				else :
					sbest=0
					count+=1
					y=dpg[i][0]
					y1=dpm[i][0]		
					df=df.append(gaiasort.loc[y,:])
					df1=df1.append(cluster1.iloc[y1,:])
					#print "gaia[",i,"][",sbest,"]=",dpg[i][sbest]
					#print "mystix[",i,"][",sbest,"]=",dpm[i][sbest]
			mini=min(magtot)
			minind=magtot.index(mini) 		
			sbest=minind
			w=dpg[i][sbest]
			w1=dpm[i][sbest]		
			df=df.append(gaiasort.loc[w,:])
			df1=df1.append(cluster1.iloc[w1,:])
			#print "gaia[",i,"][",sbest,"]=",dpg[i][sbest] 	
			#print "mystix[",i,"][",sbest,"]=",dpm[i][sbest]
		else:
			sbest=sml.index(smax)	#finding index of maximum score element in array sml
			v=dpg[i][sbest]
			v1=dpm[i][sbest]		
			df=df.append(gaiasort.loc[v,:])
			df1=df1.append(cluster1.iloc[v1,:])
			#print "gaia[",i,"][",sbest,"]=",dpg[i][sbest] 	#sbest doesn't matter
			#print "mystix[",i,"][",sbest,"]=",dpm[i][sbest]

#df=df.loc[df["astrometric_excess_noise"]<1]
df=df.reset_index(drop=True)
df1=df1.reset_index(drop=True)
result=pd.concat([df,df1],axis=1)
exp_csv=result.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/combined_data.csv',index=None,header=True)

exp_csv=df.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_gaia_data.csv',index=None,header=True)
exp_csv1=df1.to_csv(r'/home/gnr/Downloads/SC/MYStIX_GAIA/new/final_mystix_data.csv',index=None,header=True)



	 
#Number of duplicates
m=[]
for i in range(len(dpg)):
	if len(dpg[i])!=0:
		m.append(len(dpg[i]))
print len(m),max(m)

#Percenatge of duplicates
s=[]
for i in range(1,100):
	s.append(i)
#print s

count1=[]
for i in range(len(s)):
	c1=0	
	for j in m:	
		if j==s[i]:
			c1+=1
	count1.append(c1)
print count1

print count

print ww
print len(cluster1),len(gaiara),len(dpg)
