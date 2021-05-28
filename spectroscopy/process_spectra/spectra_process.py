
import os,sys,shutil,string,re
import fitsspectra as fs
import matplotlib.pyplot as plt
import getopt
import glob
import random
import subprocess

from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from scipy import integrate
from pyraf import iraf
from iraf import stsdas
from iraf import hst_calib
from iraf import synphot  
from iraf import noao  
from iraf import onedspec 
from numpy import polyfit,polyval,argmin,sqrt,mean,array,std,median
from scipy import *
from matplotlib import use as muse
from optparse import OptionParser, OptionGroup
#muse('TKagg')

#This code will run through scaling, dereddening, and deredshifting all together.
#Command line arguments: 
#	stages: -s [scaling, deredding, deredshifting, all]
#	host extinction: -he
#	galactic extinction: -ge
#	redshift: -z
#	object name: -n
#	photometricfile: -pf
#	filters: -f [UBgVrRiIzJHK]

#scaling requirements: photometric data

#1 find fits files in directory
#2 convert to ascii (need to account for all file types: based off fits header: INSTRUME)
#  
#3 

#MIGHT NEED TO APPEND TO THESE
global utKeys, mjdKeys, instKeys, parametri1
utKeys = ["DATE","DATE-OBS","DATE_OBS","DATEOBS","UTDATE"]
mdjKeys = ["ACQTIME", "JD", "MJD-OBS", "MJD", "MJD_OBS"]
instKeys = ["INSTRUME"]
parametri1={}
parametri1['ident']=['bandw','fwhm','avgwv','eqwith','zpoint']
parametri1['U']=[205.79,484.6,3652,542.62,4.327e-9] #landolt & johnson   buser and kurucz 78
parametri1['B']=[352.94,831.11,4448,1010.3,6.09e-9] #landolt & johnson   buser and kurucz 78
parametri1['V']=[351.01,826.57,5505,870.65,3.53e-9] #landolt & johnson   buser and kurucz 78
parametri1['R']=[589.71,1388.7,6555,1452.2,2.104e-9] #landolt & cousin   bessel 83
parametri1['I']=[381.67,898.77,7900.4,1226,1.157e-9] #landolt & cousin   bessel 83
parametri1['J']=[747.1,1759.3,12370,2034,3.05e-10]  #bessel in bessel and brett 88
parametri1['H']=[866.55,2040.6,16471,2882.8,1.11e-10] #bessel in bessel and brett 88
parametri1['K']=[1188.9,2799.6,22126,3664.3,3.83e-11] #bessel in bessel and brett 88
parametri1['g']= [394.17,928.19,4718.9,1158.4,5.410e-9]     # sloan (from iraf)
parametri1['r']= [344.67,811.65,6185.2,1111.2,2.495e-9]     # sloan (from iraf)
parametri1['i']= [379.57,893.82,7499.8,1044.5,1.385e-9]     # sloan (from iraf)
parametri1['z']= [514.59,1211.8,9045.6,1350.6,8.32e-10]     # sloan (from iraf)


class spectrainfo:
	def __init__(self, filename, lam, flux, hdr):
		utdate, mjd, inst = self.infofromhdr(hdr)
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

	def infofromhdr(self, hdr):
		keys = hdr.keys()
		utdate, mjd, inst = None, None, None
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
		return utdate, str(mjd), inst

class photometricexposureinfo:
	def __init__(self, f, m, m_e, mjd, n):
		self.filter = f
		self.mag = m
		self.mag_err = m_e
		self.mjd = mjd
		self.night = n

def getphotometricdata(photofile):
	lcFile_o = open(photofile).readlines()[1:]
	data = []
	for row in lcFile_o:
		rs = row.split()
		#This will only work with my scenario right now, so IDONTGIVEAFUCK
		#data.append(photometricexposureinfo(rs[3], float(rs[5]),float(rs[6]),float(rs[9]),rs[0]))
		data.append(photometricexposureinfo(rs[5], float(rs[2]),float(rs[3]),float(rs[1]),rs[0]))
	return data

def getinterpfunc_byfilter(data, filters):
	#Didn't think this would work
	#I'm a fucking genius
	interpdict = {}
	for f in filters:
		mag = [x.mag for x in data if x.filter == f]
		mjd = [x.mjd for x in data if x.filter == f]
		interp = interp1d(mjd, mag, fill_value='extrapolate')
		interpdict[f] = interp
	return interpdict

#functions from stefano's scale spectra script:
#####################
def outputdata(_filename, xdata, ydata):
    of = open(_filename, "w")
    for i,x in enumerate(xdata):
      #print "x",x, "i", i
      of.write(str(x)+" "+str(ydata[i])+"\n")
    of.close()
def residualspes(p,y,x):
	global phmin,phmax,_xerr
	err = y
	err = (y-pval(x, p,phmin,phmax))/_xerr
	#print err
	return err
def pval(x,p,pmin,pmax):
	y = 0
	for i in range(len(p)):
		y = y+p[i]*x**i
	for i in range(len(x)-1):
		if x[i]<=pmin: y[i] = y[i+1]
		if x[i]>=pmax: y[i] = y[i-1]
	return y
def DopplerMag(l,fl,b,filename,dire=None):
	#set a lower limit at 50% of the band  transmission
	#set a upper limit at 10% of the band  transmission
	bandmin = {'U':3400, 'B':4100, 'V':5000, 'R':5650, 'I':7260, 'J':10800, 'H':14500, 'K':20000,\
               'g':4200, 'r':5600, 'i':7200,'z':8000}
	bandmax = {'U':3900, 'B':5075, 'V':6520, 'R':7150, 'I':9000,  'J':13000, 'H':19000, 'K':23000,\
               'g':6520,'r':7160,'i':9000,'z':11000}


	if dire is None:
		dire='/home/samuel/Dropbox/script/scalespectra/filters/'
		#dire='/home/swyatt/Dropbox/script/scalespectra/filters/'

	#print dire
	passband={'U':dire+'landolt_u_iraf.tab',
              'B':dire+'landolt_b_iraf.tab',
              'V':dire+'landolt_v_iraf.tab',
              'R':dire+'landolt_r_iraf.tab',
              'I':dire+'landolt_i_iraf.tab',
              'g':dire+'sdss_g_005_syn.tab',
              'r':dire+'sdss_r_005_syn.tab',
              'i':dire+'sdss_i_005_syn.tab',
              'z':dire+'sdss_z_005_syn.tab',
	      'J':dire+'J_filter_persson.tab',
	      'H':dire+'H_filter_persson.tab',
	      'K':dire+'K_filter_persson.tab'}

	if l[0]>bandmin[b] or l[-1]<bandmax[b]:
		#print l[0],l[-1],bandmin[b],bandmax[b]
		return 99
	#FILE NEEDS TO BE ASCII
	#print(filename)
	#print(passband[b])
	bbb = iraf.calcphot(passband[b],spectrum=filename,form="vegamag",Stdout=1)[-1]
	#print(bbb)
	aa=string.split(bbb)[1]
	return aa
def onkeypress(event):
	import matplotlib.pyplot as plt
	from numpy import polyfit,polyval,argmin,sqrt,mean,array,std,median
	global _col,_dmag,grado,_lam_sp,_fl_sp,lines,ax,ax2,_phot_mag,_zeropoint,ppp,phmin,phmax,_xerr,fl2

	xdata,ydata = event.xdata,event.ydata
	try:
		dist = sqrt((xdata-_col)**2+(ydata-_dmag)**2)
	except:
		dist = 0
	ii = argmin(dist)
	if event.key == 'd' :
		__col,__dmag = _col.tolist(),_dmag.tolist()
		__xerr=_xerr.tolist()
		ax2.plot(_col[ii],_dmag[ii],'xk',ms=25)
		del __col[ii],__dmag[ii],__xerr[ii]
		_col,_dmag = array(__col),array(__dmag)
		_xerr=array(__xerr)

	if event.key == 'a' :
		__col,__dmag = _col.tolist(),_dmag.tolist()
		__xerr=_xerr.tolist()
		__col.append(xdata)
		__dmag.append(ydata)
		__xerr.append(1)
		_col,_dmag = array(__col),array(__dmag)
		_xerr=array(__xerr)
		ax2.plot(xdata,ydata,'db',ms=10)

	if event.key in ['1','2','3','4','5','6'] :
		grado=int(event.key)

	p0 = []
	for g in range(grado):
		p0.append(0.)
	if len(p0)>=2 and grado <= len(_dmag):
		plsq  = leastsq(residualspes,p0,args=(array(_dmag),array(_col)))
		pfit  = array(_lam_sp) # arange(phmin,phmax+1,1)
		kkfit = pval(pfit,plsq[0],phmin,phmax)
		corval = pval(_col,plsq[0],phmin,phmax)
	else:
		pfit  = array(_lam_sp) # arange(phmin,phmax+1,1)
		media=mean(_dmag)
		kkfit = array(_lam_sp)-array(_lam_sp)+media
		plsq=0

	scal=10**(array(kkfit)/2.5)
	fl2=list(array(_fl_sp)/scal)
	phot_flux=10**(array(_phot_mag)/-2.5)*array(_zeropoint)

	lines.pop(0).remove()
	lines = ax.plot(_lam_sp, fl2,'-g',label='scaled spectrum')

	ppp.pop(0).remove()
	ppp=ax2.plot(pfit,kkfit)

	plt.draw()
def splin4(ph11,fl11,fl11err,phmi,phma,lam_sp,fl_sp,zeropoint,phot_mag):
	import matplotlib.pyplot as plt
	global _col,_dmag,grado,_lam_sp,_fl_sp,lines,ax,ax2,_phot_mag,_zeropoint,ppp,_xerr,phmin,phmax,fl2

	plt.ion()
	fig = plt.figure(frameon=True)
	ax = fig.add_subplot(2,1,1,frame_on=True)
	ax2 = fig.add_subplot(2,1,2,frame_on=True)

	ph1,fl1,fl1err=[],[],[]
	for i in range(len(ph11)):
		#print phma, ph11[i], phmi
		if phma>=ph11[i]>=phmi:
			ph1.append(ph11[i])
			fl1.append(fl11[i])
			fl1err.append(fl11err[i])

	_fl_sp=fl_sp
	_lam_sp=lam_sp
	grado=1
	_col=array(ph1)
	_dmag=array(fl1)
	_xerr=array(fl1err)
	_phot_mag=phot_mag
	_zeropoint=zeropoint
	phmin=phmi
	phmax=phma

	mm=min(array(fl1))-1
	mma=max(array(fl1))+1
	phmin=phmi#(int((min(ph1)-5)*10))/10
	phmax=phma#(int(max(ph1)/10)+1)*10
	phr=phmax-phmin
	if len(str(phr))==3:
		phr=int((phmax-phmin)/100)*15
	else:
		phr=int((phmax-phmin)/10)
	if phr==0:
		phr=(phmax-phmin)/10
	mmr=abs(int((mma-mm)/10))
	if mmr==0:
		mmr=abs((mma-mm)/10)

	p0 = []
	for g in range(grado):
		p0.append(0.)

	if len(p0)>=2:
		plsq  = leastsq(residualspes,p0,args=(array(fl1),array(ph1)))
		pfit  = array(_lam_sp) # arange(phmin,phmax+1,1)
		kkfit = pval(pfit,plsq[0],phmin,phmax)
		corval = pval(_col,plsq[0],phmin,phmax)
	else:
		pfit  = array(_lam_sp) # arange(phmin,phmax+1,1)
		media=mean(fl1)
		kkfit = array(_lam_sp)-array(_lam_sp)+media
		plsq=0

	scal=10**(array(kkfit)/2.5)
	fl2=list(array(_fl_sp)/scal)
	phot_flux=10**(array(phot_mag)/-2.5)*array(_zeropoint)
	massimo=integrate.trapz(fl_sp,_lam_sp)/(max(_lam_sp)-min(_lam_sp))*4
	ax.plot(_lam_sp, _fl_sp,'-r', label='spectrum')
	lines=ax.plot(_lam_sp, fl2,'-g',label='scaled spectrum')
	lin2=ax.plot(ph11, phot_flux,'Db', label='photometry')
	ax.legend(numpoints=1,markerscale=.5)
	an = 'n'
	pp = ax2.plot(ph1,fl1,'pr')
	ppp=ax2.plot(pfit,kkfit)
	kid = fig.canvas.mpl_connect('key_press_event',onkeypress)
	#cid = fig.canvas.mpl_connect('button_press_event',onclick)
	plt.draw()
	raw_input('[a] add points, [d] remove points, [1,2,3,4] degree of poly, Return to exit ...')
	plt.close()
	return pfit,kkfit,plsq,grado,lam_sp,fl2,fl1err
#####################

#functions from stefano's deredden scripts:
####################################
def red1(x,Rv):
	####################### cardelli law foe reddeing correction
	## giving a frequency and Rv compute the absorbtion at the given frequency
	if 0.3 <= x <= 1.1:
		a=0.574*x**1.61
		b=-0.527*x**1.61
	elif 1.1 < x <= 3.3:
		y=x-1.82
		a=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
		b=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
	elif 3.3 < x <= 8:
		if x < 5.9:
			fa=fb=0
		else:
			fa=-0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
			fb=0.2130*(x-5.9)**2+0.1207*(x-5.9)**3

		a=1.752-0.316*x-0.104/((x-4.67)**2+0.341)+fa
		b=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)+fb

	Al_AV=(a+b/Rv)
	return Al_AV
def ccm(wav, rv = 3.1):
	"""
	giving one wavlength and rv gives the value to correct for reddening 
	"""
	x = 10000.0 / wav
	#print wav
	#print x
	if x < 0.3:
		print "error: wavelength is out of the range of the extinction function (1)"
	elif x < 1.1:
		y = x ** 1.61
		a = 0.574 * y
		b = -0.527 * y
	elif x < 3.3:
		y = x - 1.82
		a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 + y *(0.72085 + y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
		b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * (-5.38434 + y* (-0.62251 + y * (5.30260 + y * (-2.09002)))))))
	elif x < 5.9:
		y = (x - 4.67) ** 2
		a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
		b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
	elif x < 8.0:
		y = (x - 4.67) ** 2
		a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
		b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
		y = x - 5.9
		a = a - 0.04473 * y**2 - 0.009779 * y**3
		b = b + 0.2130 * y**2 + 0.1207 * y**3
	elif x <= 10.0:
		y = x - 8
		a = -1.072 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
		b = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3
	else:
		print "error: wavelength is out of the range of the extinctionfunction (2)"
	#rv = 3.1
	y = a + b / rv
	return y
def deredden(w_arr,d_arr,av,rv=3.1):
	"""
	if takes as input the wavelength and flux froma  spectrum 
	and correct the spectrum for the reddening
	"""
	w_arr = array(w_arr)
	d_arr = array(d_arr)
	lngth = len(w_arr)
	corr=[]
	for i in range(0, lngth):
		awavelength = w_arr[i]
		corr.append(10 ** (av * (-0.4) * ccm(awavelength,rv)))
	d_arr_dered = d_arr* corr
	return w_arr,d_arr_dered
#############################################

def main():
	description = ""
	usage = ""
	parser = OptionParser(usage=usage, description=description, version="%prog 0.5")
	parser.add_option("-a", "--runall", dest="runall", action="store_true", default=True, help="Run all fits files in the directory")
	parser.add_option("-o", "--run1", dest="runone", type=str, default=None, help="Pass only one fits file to calibrate")
	parser.add_option("-f", "--filters", dest="filters", type=str, default="", help="Filters to flux calibrate the spectra [UBgVrRiIzJHK]")
	parser.add_option("--pf", dest="photometricfile", type=str, default="", help="Ascii file containing photometric data. Assumes it is in this directory. This will probably break")
	parser.add_option("-n", "--name", dest="objectname", type=str, default="Object", help="Name of the object in the spectra")
	parser.add_option("--he", dest="ebvtoth", type=str, default=None, help="Host extinction value used to deredden spectra")
	parser.add_option("--ge", dest="ebvtotg", type=str, default=None, help="Galactic extinction value used to deredden spectra")
	parser.add_option("-z", dest="redshift", type=str, default=None, help="redshift value used to correct the spectra")
	option, args = parser.parse_args()

	filters = option.filters
	photometricdatafile = option.photometricfile
	ebvtotg = option.ebvtotg
	redshift = option.redshift
	ebvtoth = option.ebvtoth
	objectname = option.objectname
	runall = option.runall
	run1 = option.runone

        fitsshit = False

        if fitsshit:
	    if run1 is not None:
		    fitsfiles = [x for x in os.listdir(os.getcwd()) if run1 in x]
	    else: 
		    fitsfiles = [x for x in os.listdir(os.getcwd()) if ".fits" in x]
            datafiles =  fitsfiles
        else:
            if run1 is not None:
                datafiles = [x for x in os.listdir(os.getcwd()) if run1 in x]
            else:
                datafiles = [x for x in os.listdir(os.getcwd()) if '.ascii' in x]

        

	if filters is not "" and photometricdatafile is not "":
		photometricdata = getphotometricdata(photometricdatafile)
		photometricinterpfunc = getinterpfunc_byfilter(photometricdata, filters)

	for file in datafiles:
		print file

		hdu = fits.open(file)
		lam, flux = fs.getdata(hdu, file)
		hdr = hdu[0].header
		info = spectrainfo(file,lam,flux,hdr)
		if filters is not "" and photometricdatafile is not "":
			asciifile = 'asciifile_tmp.ascii'
			outputdata(asciifile, info.lam, info.flux)
			phmi,phma=min(lam),max(lam)
			dm=[]
			ref_wavelengh=[]
			zeropoint=[]
			mphot = []
			for filt in filters:
				fluxcalibdata = photometricinterpfunc[filt](info.mjd)
				mphot.append(fluxcalibdata)
				mag1=DopplerMag(lam,flux,filt,asciifile, None)
				if mag1 == 99:
					print("Selected filters outside of spectra wavelength range!")
				dm.append(fluxcalibdata-float(mag1))
				ref_wavelengh.append(parametri1[filt][2])
				zeropoint.append(parametri1[filt][-1])

			subprocess.call("rm "+asciifile, shell=True)
			zero=list(array(dm)-array(dm)+1)
		 
			phtmp,magtmp,plsq,grado,lam,fl2,fl1err=splin4(ref_wavelengh,dm,zero,phmi,phma,info.lam,info.flux,zeropoint,mphot)
			info.lam = lam
			info.flux = fl2

		if ebvtotg is not None:
			print("Correcting for galactic reddening of " + str(ebvtotg) + "\n")
			ebvtotg = float(ebvtotg)
			av=(-1)*ebvtotg*3.1
			info.lam,info.flux = deredden(info.lam,info.flux, av, 3.1)

		if redshift is not None:
			print("Correcting for redshift of " + str(redshift) + "\n")
			redshift = float(redshift)
			info.lam = array(info.lam)/(1+float(redshift))
			info.flux = array(info.flux)*(1+float(redshift))

		if ebvtoth is not None:
			print("Correcting for host reddening of " + str(ebvtotg) + "\n")
			ebvtoth = float(ebvtoth)
			av=(-1)*ebvtoth*3.1 # 3.1 
			info.lam,info.flux=deredden(info.lam,info.flux, av, 3.1)  # 3.1

		newfname = objectname + "_" + info.utdate + "_" + info.mjd + "_" + re.sub(" ", "-", info.inst) + ".ascii"
		of = open(newfname, 'w')
		for i in range(0,len(info.lam)):
			of.write('%16s\t%16s\n'  % (str(info.lam[i]),str(info.flux[i])))
		of.close()
main()
