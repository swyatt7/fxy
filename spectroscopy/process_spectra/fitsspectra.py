from astropy.io import fits
from astropy.time import Time
import os
import numpy as np
import math
import sys

def getdata(hdu, filename=""):
	hdr = hdu[0].header
	keys = hdr.keys()
	instKeys = ["INSTRUME"]
	inst = None
	for key in instKeys:
	    if key in keys:
	        inst = hdr[key]
	        if "SpeX" in inst:
	        	inst = "IRTF"
	if inst is None:
	    inst = "NOT ALFOSC"

	if "DBSP" in filename or "IMACS" in filename:
		inst = "ACAM"
	
	if "Goodman Spectro" in inst:
		inst = "SOAR"

	fl, lam = [], []
	if inst == "SOFI" or inst == "EFOSC":
		lam, fl = dataFromPESSTO(hdu,filename)
	if inst == "FIRE":
		lam, fl = dataFromFire(hdu)
	if inst == "IRTF":
		lam, fl = dataFromIRTF(hdu)
	if inst == "NOT ALFOSC":
		lam, fl = dataFromNotALFOSC(hdu)
	if inst in ['RSS', 'ACAM', 'HFOSC', 'DIS', "Bok B&C Spectrograph"]:
		lam, fl = dataFromACAM(hdu,filename)
	if inst in ['en06', 'en12', "SOAR"]:
		lam, fl = dataFromLCO(hdu, filename)
	if inst in ['Binospec']:
		lam, fl = dataFromBinospec(hdu, filename)

	return lam, fl

def dataFromPESSTO(hdu,filename):
	if	 "PSN" in filename:
		hdr = hdu[0].header
		ydata = hdu[0].data[0][0]
		naxis1 = hdr['naxis1']
		crpix1 = hdr['crpix1']
		crval1 = hdr['crval1']
		try:
			cdelt1 = hdr['cdelt1']
		except:
			cdelt1 = hdr['cd1_1']
		pix = np.array(range(1,naxis1+1,1))
		#pix = np.array(range(1,len(ydata)+1,1))
		xdata = (pix-crpix1)*cdelt1+crval1
	else:
		xdata = hdu[1].data[0][0]
		ydata = hdu[1].data[0][1]
	return xdata,ydata

def dataFromBinospec(hdu, filename):
	bunit = 1E-17
	hdr = hdu[0].header
	ydata_tmp = hdu[0].data[0]
	ydata = []
	for yy in ydata_tmp:
		ydata.append(float(yy)*bunit)
	naxis1 = hdr['naxis1']
	crpix1 = hdr['crpix1']
	crval1 = hdr['crval1']
	try:
		cdelt1 = hdr['cdelt1']
	except:
		cdelt1 = hdr['cd1_1']
	#pix = np.array(range(1,naxis1+1,1))
	pix = np.array(range(1,len(ydata)+1,1))
	xdata = (pix-crpix1)*cdelt1+crval1
	return xdata,ydata

def dataFromLCO(hdu, filename):
	hdr = hdu[0].header
	ydata = hdu[0].data[0][0]
	naxis1 = hdr['naxis1']
	crpix1 = hdr['crpix1']
	crval1 = hdr['crval1']
	try:
		cdelt1 = hdr['cdelt1']
	except:
		cdelt1 = hdr['cd1_1']
	#pix = np.array(range(1,naxis1+1,1))
	pix = np.array(range(1,len(ydata)+1,1))
	xdata = (pix-crpix1)*cdelt1+crval1
	return xdata,ydata

def dataFromACAM(hdu,filename):
	hdr = hdu[0].header
	ydata = hdu[0].data
	naxis1 = hdr['naxis1']
	crpix1 = hdr['crpix1']
	crval1 = hdr['crval1']
	try:
		cdelt1 = hdr['cdelt1']
	except:
		cdelt1 = hdr['cd1_1']
	#pix = np.array(range(1,naxis1+1,1))
	pix = np.array(range(1,len(ydata)+1,1))
	xdata = (pix-crpix1)*cdelt1+crval1
	return xdata,ydata

def dataFromFire(hdu):
    xdata = 10000*hdu[0].data[0]
    ydata = hdu[0].data[1]
    return xdata,ydata

def dataFromIRTF(hdu):
	xdata = 10000*hdu[0].data[0]
	ydata = hdu[0].data[1]
	return xdata,ydata

def dataFromNotALFOSC(hdu):
    keys = [x for x in hdu[0].header.keys() if "WAT2_" in x]
    waves = ""
    for key in keys:
        waves += hdu[0].header[key]
    #Have fun reading this shit NOTALFOSC fuckers
    swaves = waves.split('"')[1].split()[13:]
    xdata = []
    for sw in swaves:
        if len(sw.split('.')) > 2:
            swsplit = sw.split('.')
            sw1= float(swsplit[0]+'.'+swsplit[1][:-4])
            sw2= float(swsplit[1][-4:]+'.'+swsplit[2])
            xdata.append(sw1)
            xdata.append(sw2)
        else:
            xdata.append(float(sw))
    ydata = hdu[0].data

    if (xdata[0]-xdata[1] > 0):
        return xdata[::-1], ydata[::-1]
    else:
        return xdata,ydata

class SpectraInfoFits:
	def __init__(self, filename, lam, flux, hdr):
		utdate, mjd, inst, exptime = self.infofromhdr(hdr)
		if float(mjd) < 2400000.5:
			mjd = str(float(mjd) + 2400000.5)
		if len(lam) != len(flux):
			if len(lam) < len(flux):
				flux = flux[0:len(lam)]
			if len(flux) < len(lam):
				lam = lam[0:len(flux)]

		self.filename = filename
		self.lam = lam
		self.flux = flux
		self.utdate = utdate.split("T")[0]
		self.mjd = mjd
		self.inst = inst
		self.exptime = exptime

	def infofromhdr(self, hdr, filename=""):
		utKeys = ["DATE","DATE-OBS","DATE_OBS","DATEOBS","UTDATE"]
		mdjKeys = ["ACQTIME", "JD", "MJD-OBS", "MJD", "MJD_OBS"]
		instKeys = ["INSTRUME"]
		exptimeKeys = ['EXPTIME', 'EXPTOT', 'ITOT']
		keys = []
		for x in hdr.keys():
			keys.append(str(x))
		utdate, mjd, inst, exptime = None, None, None, None

		for key in exptimeKeys:
		    if key in keys:
		        exptime = hdr[key]

		for key in utKeys:
		    if key in keys:
		        utdate = hdr[key]

		for key in mdjKeys:
			if key in keys:
				if key == "ACQTIME":
					mjd = str(float(hdr[key])-2400000.5)
				else:
					mjd = hdr[key]

		for key in instKeys:
		    if key in keys:
		        inst = hdr[key]
		        if "SpeX" in inst:
		        	inst = "IRTF"

		if inst is None:
		    inst = "NOT ALFOSC"
		
		if inst == 'Goodman Spectro':
			inst = 'SOAR'

		if "DBSP" in filename:
			inst = "DBSP"
		if "IMACS" in filename:
			inst = "IMACS"

		if mjd is None and utdate is not None:
			t = Time([utdate])
			mjd = t.mjd[0]

		return utdate, str(mjd), inst, exptime
