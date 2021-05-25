import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('font', family='serif')


class exposureInfo:
	def __init__(self, f, m, m_e, mjd, n):
		self.filter = f
		self.mag = m
		self.mag_err = m_e
		self.mjd = mjd
		self.night = n

lcFile = "CSP15aak_phot_csp.dat"
lcFile_o = open(lcFile).readlines()[1:]

data = []
for row in lcFile_o:
	rs = row.split()
	data.append(exposureInfo(rs[3], float(rs[5]),float(rs[6]),float(rs[9]),rs[0]))

leg = []

offset_dict = {"u":-4.0, "B":-2.0, "g":-1.0, "V":0.0, "r":1.0, "i":1.5, "Y":3.0, "J":3.5, "H": 4.5}
filters = ['u','B','g','V','r','i','Y','J','H']
xmin = 0


#base offset for max Bmag
bdat = [x for x in data if x.filter == 'B']
peakB_i = [x.mag for x in bdat].index(min([x.mag for x in bdat]))
#2457112.71
universal_offset = 57112.71 #[x.mjd for x in bdat][peakB_i]
dist = 32.150 #taken from :http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=NGC+5839
d_mod = 5*math.log10(dist*(math.pow(10, 6))/10)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.plot((0, 0), (23, 9), "--k", alpha=0.5)

opticalspectradates = [57281.02,
					   57280.98,
					   57248.9,
					   57154.03,
					   57131.22,
                                           57136.3,
                                           57162.4,
                                           57222.0,
					   57125.11,
					   57122.34,
					   57109.38,
					   57108.16,
					   57101.32,
					   57101.29]

nirspectradates = [57135.23,
				   57149.36,
				   57128.41,
				   57174.21,
				   57159.37,
				   57131.27,
				   57131.38,
	 			   57124.22,
				   57119.29,
			   	   57114.36,
				   57110.31,
				   57102.27]
for d in opticalspectradates:
	d = d - universal_offset# + 2400000.5
	ax1.plot((d, d), (22.2,21.2), "-b")
for d in nirspectradates:
	d = d - universal_offset# + 2400000.5
	ax1.plot((d, d), (22.4,21.4), "-r")

colors = iter(cm.rainbow(np.linspace(0, 1, len(filters))))


for i,filt in enumerate(filters):
	filt_data = [x for x in data if x.filter == filt]
	mag_dat = [x.mag for x in filt_data]
	abmag_dat = [x.mag - d_mod for x in filt_data]
	mjd_dat = [x.mjd -2400000.5 for x in filt_data]

	#maxlight_i = mag_dat.index(min(mag_dat))
	#maxlightDay = mjd_dat[maxlight_i]
	mjd_datNew = [x-universal_offset for x in mjd_dat]

	offset = offset_dict[filt]
	mag_datNew = [x + offset for x in mag_dat]

	xmin = min(mjd_datNew) if min(mjd_datNew) < xmin else xmin

	strfset = str(offset) if offset < 0.0 else "+"+str(offset)
	sizes = [10 for x in mag_datNew]
	markerstyles = ['o' for x in mag_datNew]
	color = next(colors)
	if filt == 'V':
		sizes[0] = 50
		markerstyles[0] = '*'

	for it in range(0, len(mjd_datNew), 1):
		ax1.scatter([mjd_datNew[it]], [mag_datNew[it]], s=[sizes[it]], c=color, marker=markerstyles[it])
	#ax1.legend(label=filt + strfset, facecolor=color, fontsize=10)
	ax1.scatter([],[], label=filt+strfset, c=color, marker=markerstyles[1], s=8)

#ax1.legend(leg)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14}



ax1.set_xlabel(r'Days since $B_{max}$', fontdict=font)
ax1.set_ylabel(r'Apparent Magnitude', fontdict=font)
#ax1.set_title(r'Lightcurve of SN2015bp', fontdict=font)
#plt.gca().invert_yaxis()
ax2.set_ylabel("Absolute Magnitude", fontdict=font)

ax1Ys = ax1.get_yticks()
ax2Ys = []
for X in ax1Ys:
    ax2Ys.append(int(X - d_mod))


ax2.set_yticks(ax1Ys)
ax2.set_ybound(ax1.get_ybound())
ax2.set_yticklabels(ax2Ys)
ax1.set_ylim(ax1.get_ylim()[::-1])
ax2.set_ylim(ax2.get_ylim()[::-1])

ax1.legend(fontsize=9)
ax1.set_xlim(xmin-4, 100)
plt.savefig("lightcurve.png", format='png', dpi=1000)
plt.show()
#recarr = np.array([(filters), (mags), (mag_errs), (mjds), (nights)], dtype=dtype)
