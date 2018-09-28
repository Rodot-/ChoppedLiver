#!/usr/bin/env python
'''
fitsheader.py

Command line interface for viewing the header information of one or more fits files
'''

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", "--keys", help="Display Header Keys", action='store_true')
parser.add_argument("-c", "--column", help="Display Header Column", type=str)
parser.add_argument("-t", "--table", help="Print Table Header (used with -c)", action='store_true')
parser.add_argument("-u", "--hdu", help="Select Fits HDU", default=0, type=int)
parser.add_argument("files", help="Fits Files to View", type=str, nargs='+')

args = parser.parse_args()

good_files = filter(os.path.isfile, args.files)
bad_files = list(set(args.files) - set(good_files))

if len(bad_files):
	print>>sys.stderr, "The Following Provided Files Cannot be Found:"
	for f in bad_files:
		print>>sys.stderr, " ", f

assert not (bool(args.keys) and bool(args.column)), "Cannot Use Both --keys and --column"
assert not (bool(args.table) and not bool(args.column)), "Cannot Output in Table Format without Fields"

if args.column:
	keys = args.column.split(',')
	if args.table:
		print " ".join(keys)


from astropy.io import fits
for filename in args.files:
	try:
		with fits.open(filename) as f:
			header = f[args.hdu].header
			header['FILENAME'] = os.path.basename(filename)
			if args.keys:
				print " ".join(header.keys())
			elif args.column:
				print " ".join([str(header[key]) for key in keys])
			else:
				print header.__repr__()
	except Exception as e:
		print>>sys.stderr, e
