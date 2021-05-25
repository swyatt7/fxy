import fitsspectra as fs
import os
import re
from astropy.io import fits


def main():
	cwd = os.getcwd()
	fitsfiles = [f for f in os.listdir(cwd) if 'fits' in f and f[0] is not '.' and '.py' not in f]#and 'AT' in f]
	for f in fitsfiles:
		print(f)
		hdu = fits.open(f)
		lam, flux = fs.getdata(hdu, f)
		print(len(lam), len(flux))
		hdr = hdu[0].header
		info = fs.SpectraInfoFits(f,lam,flux,hdr)
		objectname = f.split('_')[0]
		newfname = objectname + "_" + info.utdate + "_" + info.mjd + "_" + re.sub(" ", "-", info.inst) + '_' + str(info.exptime) + ".ascii"
		of = open(newfname, 'w')
		for i in range(0,len(lam)):
			of.write('%16s\t%16s\n'  % (str(lam[i]),str(flux[i])))
		of.close()

main()
