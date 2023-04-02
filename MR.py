import numpy as np
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc,rcParams
from numpy import random
from astropy import constants as const

G    = const.G.value
Msun = const.M_sun.value
Rsun = const.R_sun.value
Mjup = const.M_jup.value
Rjup = const.R_jup.value
AU   = const.au.value
Mearth = const.M_earth.value
Rearth = const.R_earth.value

ecc0 = 0.42
per0 = 11.6336
mass0 = 0.18


def trans_signal(Mp,Rp,Rs,Teq,mJ):
	#https://arxiv.org/pdf/1805.03671.pdf

    try:
        n = len(Mp)
        sf = np.ones(len(Mp))
        I = np.where(Rp < 1.5)[0]
        sf[I] = 0.19
        I = np.where((Rp >=1.5) & (Rp < 2.75))[0]
        sf[I] = 1.26
        I = np.where((Rp >=2.75) & (Rp < 4))[0]
        sf[I] = 1.28
        I = np.where(Rp >=4)[0]
        sf[I] = 1.15
    except:
    	if Rp < 1.5:
    		sf = 0.19
    	elif Rp<2.75:
    		sf = 1.26
    	elif Rp < 4.:
    		sf = 1.28
    	else:
    		sf = 1.15

    return sf * Rp**3 * Teq * 10**(-mJ/5.) /(Mp * Rs**2)


# Set seaborn contexts:
sns.set_context("talk")
sns.set_style("ticks")

# Fonts:
# Arial font (pretty, not latex-like)
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
# Latex fonts, quick:
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
# Latex fonts, slow (but accurate):
rc('font', **{'family': 'Helvetica'})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size':20})
plt.rc('legend', **{'fontsize':7})

# Ticks to the outside:
rcParams['axes.linewidth'] = 3.0
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'


Qp = 3e4
Qpp = 3*Qp/(2.*0.15)

#Here goes parameters for your planet: teff,Ms,Rs,P,ecc,sma,Mp,Mp_l,Mp_u,Rp,Rp_l,Rp_u,Jmag
nd = np.array([5696.,0.95,1.123,4.375,0.,0.05146,0.09,0.014,0.014,0.614,0.020,0.020,11.184])
#Ms P ecc sma Mp Rp
#d = np.loadtxt('tepcat.txt',skiprows=0,usecols=[1,7,10,19,20,23,26,29])
d = np.loadtxt('tepcat2.txt',skiprows=0,usecols=[1,7,10,19,20,23,26,27,28,29,30,31,43])
d2 = np.loadtxt('tepcat2.txt',skiprows=0,usecols=[0],dtype='str')

d = np.vstack((d,nd))
d2 = np.hstack((d2,np.array(['NewPlanet'])))
print d
#print gfd
Mp = d[:,6]
Mpl = d[:,7]
Mpu = d[:,8]
Mpe = 0.5*(Mpl+Mpu)
Rp = d[:,9]
Rpl = d[:,10]
Rpu = d[:,11]
Rpe = 0.5*(Rpl+Rpu)
Ms = d[:,1]
sma = d[:,5]
per = d[:,3]
ecc = d[:,4]
teff = d[:,0]
Rs = d[:,2]
mJ = d[:,-1]
print Mpe/Mp

J = np.where((Mp!=-1)&(Rp!=-1)&(Mpe/Mp<0.25)&(Rp!=-1)&(Rpe/Rp<0.25)&(mJ>0))


Mp,Rp,Ms,sma,per,ecc,teff,Rs,names,mJ = Mp[J],Rp[J],Ms[J],sma[J],per[J],ecc[J],teff[J],Rs[J],d2[J],mJ[J]
print names
#print tfrd


#I = np.where((Mp>0.054)&(Mp<0.18)&(Rp>0.))
#print names[I],fds
#print I,fds
#I = np.where((ecc>0.6) & (per<30))[0]
#print I,fd

"""
A = (63./4.)*np.sqrt(G*(Ms*Msun)**3)*((Rp*Rjup)**5/(Qp*Mp*Mjup))
tau_circ = A * (sma*AU)**(-6.5)
tau_circ = 1. / tau_circ
tau_circ = tau_circ/(3600.*24*365*1e9)
I = np.where(tau_circ>100)[0]
tau_circ[I] = 100.
I = np.where(tau_circ<0.001)[0]
tau_circ[I] = 0.001

print tau_circ
#print ds
print per[-1],ecc[-1],np.log10(Mp[-1]),np.log10(tau_circ)[-1]
print np.log10(Mp)
#print fds

"""
tq = teff * np.sqrt((Rs*Rsun)/(2*sma*AU))

"""
aroche = 2.7 * Rjup*Rp[-1] * (Msun*Ms[-1]/(Mp[-1]*Mjup))**(1./3.)
tsmas = np.arange(0.01,3.,0.001)*AU
teccs = 1 - (aroche / tsmas)
tpers = np.sqrt(4.*np.pi**2*tsmas**3/(G*(Msun*Ms[-1]+Mjup*Mp[-1])))/(3600.*24.)

aroche2 = 2.7 * 1.*Rjup * (1.*Msun/(1*Mjup))**(1./3.)
tsmas2 = np.arange(0.01,3.,0.001)*AU
teccs2 = 1 - (aroche2 / tsmas2)
tpers2 = np.sqrt(4.*np.pi**2*tsmas2**3/(G*(1.*Msun+1*Mjup)))/(3600.*24.)

aroche3 = 2.7 * 0.5*Rjup * (1.*Msun/(0.05*Mjup))**(1./3.)
tsmas3 = np.arange(0.01,3.,0.001)*AU
teccs3 = 1 - (aroche3 / tsmas3)
tpers3 = np.sqrt(4.*np.pi**2*tsmas3**3/(G*(1.*Msun+0.05*Mjup)))/(3600.*24.)

H0 = (63./4.)*(sma*AU)**(-15./2.)*ecc**2*((G*Ms*Msun)**(3./2.)*Ms*Msun*(Rp*Rjup)**5)/Qpp

d3 = np.loadtxt('rv-planets.txt')

fig = plt.figure(figsize=(10,5))

plt.plot(tpers,teccs,'r',alpha=0.3)
plt.plot(tpers2,teccs2,'k',alpha=0.7)
plt.plot(tpers3,teccs3,'b',alpha=0.7)

plt.scatter(d3[:,1],d3[:,2],s=25*np.log10(d3[:,0]*3178.),c='black',alpha=0.1)
plt.scatter(per,ecc,s=25*np.log10(Mp*3178.),c=tau_circ,norm=matplotlib.colors.LogNorm())

plt.plot([nd[3],nd[3]],[nd[4]+0.03,nd[4]+0.1],'g',alpha=0.6)
plt.plot([nd[3],nd[3]],[nd[4]-0.03,nd[4]-0.1],'g',alpha=0.6)
plt.colorbar(label=r'$\tau_{circ}$ [Gyr]')

x=[0.5,0.5,0.5,0.5]
y=[0.7,0.76,0.82,0.88]
s=[10.,170.,3180,3180*10]

plt.scatter(x,y,s=25*np.log10(s),c='black')
plt.text(0.6,0.68,r'1 M$_{\oplus}$', fontsize=15)
plt.text(0.6,0.74,r'1 M$_{N}$', fontsize=15)
plt.text(0.6,0.8,r'1 M$_{J}$', fontsize=15)
plt.text(0.6,0.86,r'10 M$_{J}$', fontsize=15)

plt.xscale("log", nonposx='clip')
plt.xlabel('Period  [d]', fontsize=18)
plt.ylabel('eccentricity', fontsize=18)
plt.xlim([0.1,2000])
plt.ylim([0.,1])
plt.tick_params(labelsize=20)

plt.tight_layout()
plt.savefig("eccs.pdf", dpi=300)

plt.clf()
"""

ts = trans_signal(Mjup*Mp/Mearth,Rjup*Rp/Rearth,Rs,tq,mJ)
#print ts

d = np.loadtxt('fortney.txt')
fig = plt.figure(figsize=(10,5))
plt.axvspan(0.054, 0.180, facecolor='g', alpha=0.1)
plt.xscale("log", nonposx='clip')
plt.scatter(Mp,Rp,s=ts/2.,c=tq,alpha=0.8)
plt.plot([Mp[-1]-Mpl[-1],Mp[-1]-0.05],[Rp[-1],Rp[-1]],'b',alpha=0.7)
plt.plot([Mp[-1]+0.05,Mp[-1]+Mpu[-1]],[Rp[-1],Rp[-1]],'b',alpha=0.7)
plt.plot(d[:,0],d[:,1],'k',alpha=0.7)
plt.colorbar(label=r'T$_{eq}$ [K]')

x=[0.02,0.02,0.02]
y=[1.53,1.62,1.76]
s=0.5*np.array([10.,100.,1000])

plt.scatter(x,y,s=s,c='black')
plt.text(0.025,1.51,r'10', fontsize=15)
plt.text(0.025,1.60,r'100', fontsize=15)
plt.text(0.025,1.74,r'1000', fontsize=15)
plt.text(0.02,1.90,r'TSM', fontsize=15)

tms   = np.arange(0.01,10,0.001)
iso1  = ((tms*Mjup)/(4*np.pi*300./3.))**(1./3)
iso2  = ((tms*Mjup)/(4*np.pi*3000/3.))**(1./3)
iso3  = ((tms*Mjup)/(4*np.pi*30000/3.))**(1./3)

plt.plot(tms,iso1/Rjup,'k--',alpha=0.2)
plt.plot(tms,iso2/Rjup,'k--',alpha=0.2)
plt.plot(tms,iso3/Rjup,'k--',alpha=0.2)

plt.xlabel(r'M$_P$ [M$_J$]')
plt.ylabel(r'R$_P$ [R$_J$]')
plt.xlim([0.01,10])
plt.ylim([0.,2.2])
plt.tick_params(labelsize=20)
plt.tight_layout()
plt.savefig("m-r.pdf", dpi=300)

#I = np.where((ecc>0.4) & (per<50))[0]
#for i in I:
#    print names[i],Mp[i],Rp[i],ts[i], ecc[i], per[i],tau_circ[i]
