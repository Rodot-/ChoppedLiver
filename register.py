#!/usr/bin/env python
'''
register.py

Determine which field belongs to which image
and produce a data file containing information
about each fits file

'''
import glob
from astropy.io import fits
import numpy as np
from scipy.cluster.vq import kmeans2
from matplotlib.pyplot import subplots, show
import os
import sys

path = '../ISIS_scripts/optical-data/fits_files/'
if len(sys.argv) > 1:
	path = sys.argv[1]

images = glob.glob(os.path.join(path, '*.fits'))

def get_coords(fname):

	with fits.open(fname) as f:

		header = f[0].header

		ra_hms = header['RA']
		dec_dms = header['DEC']
		ra_coords = map(float, ra_hms.split(':'))	
		dec_coords = map(float, dec_dms.split(':'))	

		ra = ra_coords[0]*15+ra_coords[1]*15/60.0+ra_coords[2]*15/3600.0
		dec = dec_coords[0]+dec_coords[1]/60.0+dec_coords[2]/3600.0

	return ra, dec

def get_header_info(fname):

	keys = ['MJD','SEEING','RA','DEC','GAIN','READNOIS','DATE-OBS']
	print(fname)
	with fits.open(fname) as f:
		return [str(f[0].header[k]) if k in f[0].header else 'None' for k in keys]

coords = np.array(map(get_coords, images))

init = np.array([[150.2, 2.5],
				[150.35, 2.5],
				[150.5, 2.5],
				[150.2, 2.35],
				[150.35, 2.35],
				[150.5, 2.35],
				[150.2, 2.2],
				[150.35, 2.2],
				[150.5, 2.2]]) # initial guesses to field centroids
				

centroids, labels = kmeans2(coords, init, iter=20, minit='matrix')

with open('all_obs.dat','w') as f:

	rows = []
	header='file,group,MJD,seeing,ra,dec,gain,readnoise,date'
	rows.append(header)
	for fname, label in zip(images, labels):

		info = get_header_info(fname)
		fname = os.path.basename(fname)
		row = ','.join([fname,str(label)]+info)
		rows.append(row)
	f.write('\n'.join(rows))

fig, ax = subplots(1,1)
for i in set(labels):
	ax.scatter(*zip(*coords[labels == i]))
ax.set_xlabel('RA')
ax.set_ylabel('DEC')
show()


