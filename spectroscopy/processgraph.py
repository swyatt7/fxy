import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import os
import numpy as np
import math
import sys
import re
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('font', family='serif')


global maxB
maxB = 0

global explos
explos = maxB - 19.08

def UTdatefromMJD(mjd):
    t = Time(mjd, format='mjd', scale='utc')
    timedec = (t.datetime.hour * 3600 + t.datetime.minute * 60 + t.datetime.second)/86400.
    strdatetime = str(t.datetime).split()[0]+'.'+str(round(timedec,2)).split('.')[1]
    if len(strdatetime) == 12:
        strdatetime = strdatetime + "0"
    return strdatetime

class spec_table:
    def __init__(self, utdate, mjd, inst, ittime):
        if mjd > 2400000:
            mjd = mjd - 2400000.5
        self.utdate = UTdatefromMJD(mjd)
        self.mjd = str(round(float(mjd), 1))
        self.inst = inst
        self.tmaxB = str(round(float(mjd) - maxB, 1))
        self.texp = str(round(float(mjd) - explos, 1))
        self.ittime = str(round(float(ittime)/60.0, 1))

    def filestring_tex_table(self):
        return "{} & {} & {} & {} & {}".format(self.utdate, self.mjd, self.inst, self.tmaxB, self.ittime)

def latex_table(sorted_dates, datedictionary, data):
    out = "\\begin{deluxetable*}{ccccc}\n"
    out += "\\tablecaption{Journal of the "+data+" spectroscopic observations \\label{tab:params}}\n"
    out += "\\tablehead{\\colhead{UT Date} & \\colhead{MJD} & \\colhead{Instrument} & \\colhead{$t_{max}(B)$} & \\colhead{$t_{exp}$} &\\colhed{$T_{int}$}\n"
    out += "\\startdata\n"
    for date in sorted_dates:
        file = datedictionary[date]
        info = spec_table_info_fromfile(file)
        out += info.filestring_tex_table() + " \\\\\n"
    out += "\\enddata\n"
    out += "\\end{deluxetable*}\n"
    of = open(data+"_tex_table.txt", "w")
    of.write(out)
    of.close()

def spec_table_info_fromfile(_file):
    _file = _file.split('PROCESS')[1]
    _file = _file.split('.p')[0]
    print(_file)
    sf = _file.split('.ascii')[0].split('_')
    objectname = sf[0]
    utdate = sf[1]
    mjd = float(sf[2])
    inst = re.sub('-', "", sf[3])
    if 'NOT' in inst:
        inst = re.sub('NOT', '', inst)
    if 'HCT' in inst:
	    inst = re.sub('HCT', '', inst)
    if 'en06' in inst:
	    inst = 'OGG'
    if 'en12' in inst:
	    inst = 'COJ'
    ittime = sf[4].split('.')[0]
    return spec_table(utdate,mjd,inst,ittime)

def getutdatefromfile(_file):
    return _file.split('_')[2]

def dataFromAsci(ascifile):                              # read ascii file
    f=open(ascifile,'r')
    s=f.readlines()
    f.close()
    vecb1,vecb2 =  [],[]
    for x in s:
        if x[0]!="#":
            try:
                c1,c2 = x.split()
            except:
                c1,c2 = x.split(',')
            vecb1.append(float(c1))
            vecb2.append(float(c2))
        if '20150412' in ascifile:
            vecb1 = vecb1[::-1]
    return vecb1,vecb2

files = [x for x in os.listdir(os.getcwd()) if ".ascii.p" in x]
datedictionary = {}
for f in files:
	date = getutdatefromfile(f)
	datedictionary[date] = f

sorted_dates_rev = sorted(datedictionary.keys(), reverse=True)
sorted_dates = sorted(datedictionary.keys())
jds = []
latex_table(sorted_dates, datedictionary, "optical")
#sys.exit()
offset = 0
xmin, xmax, ymin, ymax = 0, 0, 0, 0
figsize = plt.figaspect(0.75/1.0)

#fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
fig = plt.figure(figsize=figsize)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
yminoverall = 0
ymaxprev = 0
offsetdict = {0:15., 1:16., 2:17, 3:15.65, 4:17.2, 5:17.77, 6:18.3, 7:18.82, 8:19.5, 9:20, 10:20.5, 11:21, 12:21.8, 13:22.5, 14:24}
annpos = {0:-1.13, 1:-0.7288, 2:0.386, 3:1.62, 4:2.32, 5:2.86, 6:3.34, 7:3.85, 8:4.333, 9:4.9, 10:5.4, 11:6.1, 12:6.7, 13:7.38, 14:7.38}

specfiles = sorted_dates_rev
avgflux = 0
flux_thresh = 0.5
yval_thresh = 0
for ii,f in enumerate(specfiles):
    avgflux = avgflux+flux_thresh

    _file = datedictionary[f]
    print(_file)
    info = spec_table_info_fromfile(_file)

    print(ii,_file)

    instrument = ""

    xdata,ydata = dataFromAsci(_file)

    newx, newy = [], []

    itera = 0
    for yval in ydata[:len(xdata)]:
        if yval > 0 and xdata[itera] < 9000 and xdata[itera] > 3600:
            newx.append(xdata[itera])
            newy.append(math.log10(yval) + offset)
        itera += 1

    if float(info.tmaxB) > 100:
        flux_thresh =  0.7
    else:
        flux_thresh = 0.5
    thisavg = np.mean(newy)
    offset = avgflux-thisavg-flux_thresh
    newy = [ny + offset for ny in newy]

    yavg = np.mean(newy)
    ystd = np.std(newy)
    nnx, nny = [],[]
    if float(info.tmaxB) > 100:
	    for iterr,yeetval in enumerate(newy):
	        if abs(yavg - yeetval) < 3.5*ystd:
	            nnx.append(newx[iterr])
	            nny.append(yeetval)
	    newx = nnx
	    newy = nny
    if xmin == 0:
        xmin = min(newx)

    ymin, ymax = min(newy), max(newy)

    xmin = min(newx) if xmin > min(newx) else xmin
    xmax = max(newx) if xmax < max(newx) else xmax

    yminoverall = min(newy) if min(newy) < yminoverall else yminoverall

    ax1.plot(newx, newy, "-k", linewidth=0.666)
    anpos = newy[len(newy)-1]+yval_thresh
    ax1.annotate(str(info.tmaxB) + " - " + info.inst, xy=(9300, anpos), fontsize = 6)
    jds.append(info.tmaxB)
print()
print(jds)
ylabel = r'$\rm{log}_{10}(\rm{F}_{\lambda})\rm{+constant}$'
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 11}

titleFont = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16}

ax1.set_xlabel(r'Observed Wavelength $(\rm{\AA})$', fontdict=font)
ax2.set_xlabel(r'Rest Wavelength $(\rm{\AA})$', fontdict=font)

ax1.set_xlim(3550, 10600)
#ax1.set_ylim(-2, 10) #nice
ax1Xs = ax1.get_xticks()
ax2Xs = []
z = 0.009393
for X in ax1Xs:
    ax2Xs.append(round(X/(z+1),1))


ax2.set_xlim(3550.0/(1+z), 10600.0/(1+z))
#ax2.set_xticklabels(ax1Xs/10000)
#ax2.set_xticks(ax1Xs)
#ax2.set_xbound(ax1.get_xbound())
#ax2.set_xticklabels(ax2Xs)
#Optical Skyline 6850-6960 - center 6905
#                7540-7700

ax1.axvline(x=6905, linewidth=5, color='grey', alpha=0.5)
ax1.axvline(x=7600, linewidth=10, color='grey', alpha=0.5)

telpos = np.mean(newy) + 0.4

ax1.scatter(6905,telpos, marker='o', c='none', edgecolor='k', s=70)
ax1.scatter(6905,telpos, marker='+', c='k', edgecolors='none', s=70)
ax1.scatter(7600,telpos, marker='o', c='none', edgecolor='k', s=70)
ax1.scatter(7600,telpos, marker='+', c='k', edgecolors='none', s=70)


#ax1.legend(loc='lower left', fontsize=8)
ax1.set_ylabel(ylabel, fontdict=font)
#title = ax1.set_title(r'Spectral Evolution of SN 2015bp', fontdict=titleFont)
#title.set_y(1.1)
fig.subplots_adjust(top=0.85)

plt.savefig("Optical-Spectral_evolution.png", format='png', dpi=200)
plt.show()
