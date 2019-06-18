'''
median stack some fits images
'''

import sys
from astropy.io import fits
import numpy as np

from matplotlib.pyplot import subplots, show

def main():

	if len(sys.argv) == 1:
		print("Please select files to stack")
	files = sys.argv[1:]	

	arrays = []
	for i,f in enumerate(files):
		with fits.open(f) as ff:
			arrays.append(np.array(ff[0].data))
	arrays = np.array(arrays)

	out = np.nanmedian(arrays, axis=0)
	#fig, ax = subplots(1,1)
	#ax.imshow(np.log10(np.abs(out)))
	var = np.nanvar(arrays, axis=0)

	out_hdu = fits.PrimaryHDU(out)
	var_hdu = fits.PrimaryHDU(var)

	fits.HDUList([out_hdu]).writeto('stacked.fits', overwrite=True)
	fits.HDUList([var_hdu]).writeto('var.fits', overwrite=True)


if __name__ == '__main__':
	
	main()

	show()
