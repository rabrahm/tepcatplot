import os
import urllib
from lxml import html
from bs4 import BeautifulSoup

import astropy.coordinates as coord
from astroquery.irsa import Irsa
import astropy.units as u
from astropy.coordinates import SkyCoord

import numpy as np

def ratodeg(sra):
	ras = sra.split()
	h = float(ras[0])*360./24.
	m = (float(ras[1])/60.)*360./24.
	s = (float(ras[2])/3600.)*360./24.
	return h+m+s
def dectodeg(sdec):
	neg = False
	decs = sdec.split()
	try:
		d = float(decs[0])
	except:
		d = -float(decs[0][-2:])
		neg = True
	m = float(decs[1])/60.
	s = float(decs[2])/3600.
	if neg:
		return d - m - s
	else:
		return d + m + s

def get_page(planet):
	base_url = 'http://www.astro.keele.ac.uk/jkt/tepcat/planets/'
	query_url = "{0}{1}.html".format(base_url, planet)
	html = urllib.urlopen(query_url).read()
	soup = BeautifulSoup(html, 'html.parser')
	return soup


os.system('wget http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-ascii.txt')
planets = np.loadtxt('allplanets-ascii.txt',skiprows=0,usecols=0,dtype='str')

lines = open('allplanets-ascii.txt').readlines()
line0 = lines[0]
lines = lines[1:]
i=0
for line in planets:
	planet = line.split()[0]
	soup = str(get_page(planet))
	ra = soup.split('<tr bgcolor="#ffffff"> <td align="LEFT" style="WHITE-SPACE: NOWRAP"> Right ascension </td> <td align="CENTER" style="WHITE-SPACE: NOWRAP">')[1][:12]
	dec= soup.split('<tr bgcolor="#ffffcc"> <td align="LEFT" style="WHITE-SPACE: NOWRAP"> Declination </td> <td align="CENTER" style="WHITE-SPACE: NOWRAP">')[1][:15]
	dec = dec.split('<')[0]
	RA,DEC = ratodeg(ra),dectodeg(dec)
	#print planet,ra, dec,RA,DEC
	table1 = Irsa.query_region(coord.SkyCoord(RA, DEC, unit=(u.deg,u.deg)), catalog="fp_psc", spatial="Cone",radius=5*u.arcsec)
	vec = np.array(table1['j_m'])
	if len(vec)>0:
		print planet, vec[0]
		lines[i] = lines[i][:-1]+'\t'+str(vec[0])+'\n'
	else:
		print planet, -1
		lines[i] = lines[i][:-1]+'\t'+str(-1)+'\n'
	print lines[i]
	i+=1
f=open('tepcat2.txt','w')
f.writelines(lines)
f.close()


