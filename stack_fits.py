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

	arrays = np.empty((len(files), 2048, 2048), dtype=np.float)
	for i,f in enumerate(files):
		with fits.open(f) as ff:
			arrays[i] = np.array(ff[0].data)

	out = np.nanmedian(arrays, axis=0)
	#fig, ax = subplots(1,1)
	#ax.imshow(np.log10(np.abs(out)))
	var = np.nanvar(arrays, axis=0)

	out_hdu = fits.PrimaryHDU(out)
	var_hdu = fits.PrimaryHDU(var)

	fits.HDUList([out_hdu]).writeto('stacked.fits')
	fits.HDUList([var_hdu]).writeto('var.fits')


if __name__ == '__main__':
	
	main()

	show()
