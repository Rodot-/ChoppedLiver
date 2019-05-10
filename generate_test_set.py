'''
generate_test_set.py

A script for cleaning and preparing 
optical images of the CHILES VERDE field
from the Liverpool Telescope for processing
with ISIS

Create a set of images (and dates files)
for processing using ISIS according to some parameter.

At the time of creation, the parameter is objects around 70 degrees rotation

Author - Jack O'Brien
Date - 7/23/2018
Modified - 8/10/2018
'''

__author__ = "Jack O\'Brien"
__version__ = "0.0.1"
__maintainer__ = "Jack O'Brien"
__email__ = "obrie278@msu.edu"
__date__ = "8/10/2018"

import csv
import sys
import os
import shutil
import glob
from datetime import date

from PIL import Image
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
#from reproject import reproject_interp

import cosmics
from context import LoadingBar

# These must be set by the user
DATA_DIR = '../ISIS_scripts/data/clean' # This is where we put our cleaned images
MASK_DIR = '../ISIS_scripts/data/mask' # This is where we put our image masks
ISIS_DIR = '../ISIS_scripts/package' # This is the ISIS installation directory
RAW_DATA_DIR = '../ISIS_scripts/optical-data/fits_files' # this is where the original fits images are


INTERPOLATION = Image.BICUBIC # slowest but best rotation interp
DATA_FILE = 'cleandatatable.dat' # this will be created
TARGET_DIR = 'optical_reduction' # This is the directory containing our ISIS directories
DATES_FILE = os.path.join(TARGET_DIR, 'register2', 'dates') # depricated
DEFAULT_KEYS = ("FILENAME", 
				"MJD", # Date
				"SEEING", # PSF Sigma
				"OBSID", # Can only select OBS_1 for now
				)

# This should be a more sophisticated function that gets 
# observation numbers from a data table
# right now we're only working on obs_1
SELECT = lambda row: row['OBSID'] == 'OBS_1'

def mkdir(directory, clear=False):

	#head = os.path.split(directory)[0]
	#if head != directory and head != '':
	#	return mkdir(head)
	if not os.path.isdir(directory):
		os.mkdir(directory)

def setup_ISIS(run_name, **kwargs):
	'''set up some ISIS files in the dir to be run in'''

	register_dir = os.path.join(TARGET_DIR, 'register_{}'.format(run_name))
	images_dir = os.path.join(TARGET_DIR, 'images_{}'.format(run_name))

	mkdir(TARGET_DIR)
	
	mkdir(register_dir)
	print "created {}".format(register_dir)
	mkdir(images_dir)
	print "created {}".format(images_dir)

	dates_file = os.path.join(register_dir, 'dates')

	rows = sorted(get_rows(), key=lambda x: float(x[2]))
	with open(dates_file, 'w') as f:
		f.write('\n'.join([' '.join(row) for row in rows]))	
	print "wrote dates file"

	print "Copying Image Files..."
	for row in LoadingBar(True).load(rows): # copy over the image files
		src = os.path.join(DATA_DIR, row[0])
		dst = os.path.join(images_dir, row[0])
		shutil.copyfile(src, dst)
		#print "copied {} to {}...".format(src, dst)


	print "Copying Scripts and Config Files..."
	scripts = glob.glob(os.path.join(ISIS_DIR,'register2/*.csh'))
	scripts.extend(glob.glob(os.path.join(ISIS_DIR,'register2/*config')))
	for src in LoadingBar(True).load(scripts): # copy over the scripts
		dst = os.path.join(register_dir, os.path.basename(src))
		shutil.copyfile(src, dst)
		#shutil.copymode(src, dst) # so they remain executable
		#print "copied {} to {}".format(src, dst)

	pc_fields = ["IM_DIR", "MRJ_DIR", "REFERENCE", "REF_SUB",
				"INFILE", "VARIABLES", "DEGREE", "CONFIG_DIR",
				"SIG_THRESH", "COSMIC_THRESH", "REF_STACK",
				"N_REJECT", "MESH_SMOOTH"]

	pc_comments = {
		"IM_DIR":"  Directory where the images are",
		"MRJ_DIR":"  Installation directory",
		"REFERENCE":"  Reference image for astrometry",
		"REF_SUB":"  Reference image for subtraction",
		"INFILE":"  Dates of the frames",
		"VARIABLES":"  coordinates of objects for which we want to make light curves",
		"DEGREE":"  The degree of the polynomial astrometric transform between frames",
		"CONFIG_DIR":"  Where to find the configuration files",
		"SIG_THRESH":" ",
		"COSMIC_THRESH":" ",
		"REF_STACK":" ",
		"N_REJECT":" ",
		"MESH_SMOOTH":" ",
		}
	pc = {}

	pc["IM_DIR"] = "../images_{}".format(run_name)
	pc["MRJ_DIR"] = "../../package"
	pc["REFERENCE"] = "{}".format(rows[0][1])
	pc["REF_SUB"] = "ref.fits"
	pc["INFILE"] = "../register_{}/dates".format(run_name)
	pc["VARIABLES"] = "phot.data"
	pc["DEGREE"] = "1" # this could be edited to kwargs maybe?
	pc["CONFIG_DIR"] = "../register_{}".format(run_name)
	pc["SIG_THRESH"] = "0.27" # ditto
	pc["COSMIC_THRESH"] = "1000.0" # ditto
	pc["REF_STACK"] = "{}".format(rows[0][1])
	pc["N_REJECT"] = "2" #ditto
	pc["MESH_SMOOTH"] = "3" #ditto

	for key in pc_fields: # allow custom parameters in kwargs
		if key in kwargs:
			pc[key] = str(kwargs[key])

	config_file = os.path.join(register_dir, 'process_config')
	with open(config_file, 'w') as f:
		for field in pc_fields:
			line = (field, pc[field], pc_comments[field], '\n')
			f.write(' '.join(line))
	print "wrote out process_config file"
	

def clean_image(infile, true_rot=True, crop=True, wcs=None):

	basefile = os.path.basename(infile)
	cleanfile = os.path.join(DATA_DIR, basefile.replace('.fits','_clean.fits'))
	maskfile = os.path.join(MASK_DIR, basefile.replace('.fits','_mask.fits'))
	# Read the FITS :
	array, header = cosmics.fromfits(infile)
	# array is a 2D numpy array
	if wcs is None:
		try:
			clean_image.wcs
		except AttributeError:
			#clean_image.wcs = WCS(header) # use the wcs coord of the first input
			clean_image.wcs = header # use the wcs coord of the first input

	# crop array to be a uniform shape (2048x2048 square)

	# check the rotation angle on the sky

	wcs = WCS(clean_image.wcs)
	rotangle = header['ROTSKYPA']

	# perform WCS transform
	if true_rot: # use PIL to rotate the image
		hdu = fits.open(infile)[0]
		array = reproject_interp(hdu, clean_image.wcs)[0] 
		print "reprojected to", WCS(clean_image.wcs)
		#im = Image.fromarray(array)
		#if header['CCDRDOUT'] == 'DUMMYBOTTOMRIGHT':
		#	rotangle = rotangle+180
		#rot_im = im.rotate(rotangle, expand=True, 
		#						resample=INTERPOLATION)
		#array = np.array(rot_im)

	else: # use numpy to rotate to nearest 90 degrees
		# find the number of requred rotations
		rotangle = int(-rotangle/10+0.5)
		k = (rotangle/9)
		# make the number positive if necessary
		k = (k+4) % 4
		# and finally rotate the array
		array = np.rot90(array, k=k, axes=(0,1))
	
	if crop: # Helps ISIS interpolate more accurately
		array = array[:2048,:2048]

	# get parameters
	gain = header['GAIN']
	readnoise = header['READNOIS']
	sigclip = 5.0
	sigfrac = 0.3
	objlim = 5.0
	maxiter=4

	# Build the object :
	c = cosmics.cosmicsimage(array, gain=gain, 
							readnoise=readnoise, sigclip=sigclip, 
							sigfrac=sigfrac, objlim=objlim)
	# There are other options, check the manual...

	# Run the full artillery :
	c.run(maxiter=maxiter)

	
	header.update(wcs.to_header())

	# update our header to describe cleaning
	header['IC_PREF'] = 'IC CP'
	header.comments['IC_PREF'] = 'IC: Image Correct Params, CP: Cosmics Params'

	header['IC_DATE'] = date.today().isoformat()
	header.comments['IC_DATE'] = 'Date corrections were made'
	header['IC_VER'] = __version__
	header.comments['IC_VER'] = 'Correction script version number'
	header['IC_SCRPT'] = os.path.basename(__file__)
	header.comments['IC_SCRPT'] = 'Name of script used for corrections'
	header['IC_ATHR'] = __author__
	header.comments['IC_ATHR'] = 'Correction script author'
	header['IC_MOD'] = __date__
	header.comments['IC_MOD'] = 'Date of last modification to correction script'
	header['IC_MNTNR'] = __maintainer__
	header.comments['IC_MNTNR'] = 'Correction script maintainer'
	header['IC_EMAIL'] = __email__
	header.comments['IC_EMAIL'] = 'Email of current correction script maintainer'
	#header['IC_ROT'] = 'PIL.Image.rotate' if true_rot else 'numpy.rot90'
	header['IC_ROT'] = 'WCS' if true_rot else 'numpy.rot90'
	header.comments['IC_ROT'] = 'Rotation, PIL: exact, numpy: nearest 90 deg'
	header['IC_CROP'] = crop
	header.comments['IC_CROP'] = 'If the Image was cropped to 2048x2048'
	header['IC_INTRP'] = 'BICUBIC' if true_rot else 'None'
	header.comments['IC_INTRP'] = 'PIL rotation interpolation method'

	header['CP_COSMC'] = "L.A.Cosmic"
	header.comments['CP_COSMC'] = 'Package used for cosmic ray removal'
	header['CP_SGCLP'] = sigclip
	header.comments['CP_SGCLP'] = 'sigclip kwarg for cosmic ray removal'
	header['CP_SGFRC'] = sigfrac
	header.comments['CP_SGFRC'] = 'sigfrac kwarg for cosmic ray removal'
	header['CP_OBJLM'] = objlim
	header.comments['CP_OBJLM'] = 'objlim kwarg for cosmic ray removal'
	header['CP_GAIN'] = gain
	header.comments['CP_GAIN'] = 'gain kwarg for cosmic ray removal'
	header['CP_RDNSE'] = readnoise
	header.comments['CP_RDNSE'] = 'readnoise kwarg for cosmic ray removal'
	header['CP_MITR'] = maxiter
	header.comments['CP_MITR'] = 'maxiter kwarg for cosmic ray removal'

	# Write the cleaned image into a new FITS file, conserving the original header :
	cosmics.tofits(cleanfile, c.cleanarray, hdr=header, cols={'MASK':c.mask})

	# (c.mask is a boolean numpy array, that gets converted here to an integer array)


def clean_images():

	fits_files = glob.glob(os.path.join(RAW_DATA_DIR, '*.fits'))
	cleanable_files = []

	###
	max_clean=np.inf
	###

	for ff in fits_files:
		with fits.open(ff) as f:
			if f[0].header['OBSID'] == 'OBS_1':
				cleanable_files.append(ff)

	for f in cleanable_files:
		clean_image(f)

def generate_data_table(keys):

	arg = ",".join(keys)
	stream = os.popen('./fitsheader.py -t -c {} {}/*.fits'.format(arg, DATA_DIR))
	with open(DATA_FILE, 'w') as f:
		f.write(stream.read())
	print "created {}".format(DATA_FILE)

def get_rows(data_file=DATA_FILE):
	with open(data_file, 'r') as df:
		table = csv.DictReader(df, delimiter=' ')
		for row in table:
			if row['OBSID'] == 'OBS_1':
				#if os.path.isfile(os.path.join(DATA_DIR, row['FILENAME'].replace('.fits','_clean.fits'))):
				yield (row['FILENAME'], row['MJD'], row['SEEING'])	

def reset(run_name):

	old_files = [DATA_FILE]
	for f in old_files:
		try:
			os.remove(f)
		except OSError as e:
			print e

	register_dir = os.path.join(TARGET_DIR, 'register_{}'.format(run_name))
	images_dir = os.path.join(TARGET_DIR, 'images_{}'.format(run_name))
	old_dirs = [register_dir, images_dir]	
	for d in old_dirs:
		try:
			shutil.rmtree(d)
		except OSError as e:
			print e
	old_files = os.popen('ls %s | grep -v "h_e*_clean.fits"' % images_dir)

def test():

	import glob
	files = glob.glob('optical-data/fits_files/*.fits')
	clean_image(files[0])

def main():

	run_name = 'obs1'
	args = sys.argv[1:]
	if '-d' in args:
		args.extend(['-c','-t','-i'])
		
	if '--reset' in args:
		reset(run_name)

	if ('-i' in args) and ('-t' not in args):
		if not os.path.isfile(DATA_FILE):
			args.append('-t')
	if '-h' in args or '--help' in args:
		print "Script for managing and setting up optical data"
		print " "
		print "-h, --help: This Message"
		print "--reset:    Reset files and directories"
		print "-c:         Clean Images"
		print "-t:         Create Data Tables"
		print "-i:         Setup ISIS directories"
		print "-d:         Same as \"-c -t -i\""
		sys.exit(0)

	if '--test' in args:
		test()
		exit(0)

	if '-c' in args:
		print "Cleaning"
		clean_images()
	if '-t' in args:
		print "Creating data table"
		generate_data_table(DEFAULT_KEYS)
	if '-i' in args:
		print "Setting up ISIS Dirs"
		setup_ISIS(run_name)

if __name__ == '__main__':

	main()
