import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('font', family='serif')


class exposureInfo:
	def __init__(self, f, m, m_e, mjd, n, fname):
		self.fname = fname
		self.filter = f
		self.mag = m
		self.mag_err = m_e
		if mjd < 60000:
			self.mjd = mjd
		else:
			self.mjd = mjd - 2400000.5
		self.night = n

class lcfile():
	def __init__(self, fname, jd_index, filter_index, mag_index, mag_err_index, night_index=None, delim=None, headerrow=False):
		self.fname = fname
		self.delim = delim
		self.jd_index = jd_index
		self.filter_index = filter_index
		self.mag_index = mag_index
		self.mag_err_index = mag_err_index
		self.night_index = night_index
		self.headerrow = headerrow
		self.data = []

	def ingest(self):
		if self.headerrow:
			lcfile_o = lcFile_o = open(self.fname).readlines()[1:]
		else:
			lcfile_o = lcFile_o = open(self.fname).readlines()
		self.data = []
		for row in lcFile_o:
			rs = row.split(self.delim)
			if '#' not in row:
				if self.night_index is not None:
					night = rs[self.night_index]
				else:
					night = ''
				
				if self.filter_index == -1:
					filt='clear'
				else:
					filt=rs[self.filter_index]

				self.data.append(
					exposureInfo(
						filt,
						float(rs[self.mag_index]),
						float(rs[self.mag_err_index]),
						float(rs[self.jd_index]),
						night, 
						self.fname 
					)
				)
		return self.data

def output_for_snpy(snname, z, ra, dec, overall_data):
	fname = '{}_snpy.txt'.format(snname)
	filters = list(set([f.filter for f in overall_data]))
	filters = ['U', 'B', 'V', 'g', 'r', 'i', 'z']
	filter_snpy = {'U':'Us', 'B':'Bs', 'V':'Vs', 'g':'g_s', 'r':'r_s', 'i':'i_s', 'z':'z_s'}
	with open(fname, 'w') as of:
		of.write('{} {} {} {}\n'.format(snname, z, ra, dec))
		for f in filters:
			filtdata = [x for x in overall_data if x.filter == f]
			if len(filtdata):
				of.write('filter {}\n'.format(filter_snpy[f]))
				for fd in filtdata:
					of.write('{} {} {}\n'.format(fd.mjd, fd.mag, fd.mag_err))


def mean_mjd_data(data, day_thresh=0.4):
	mmd = []
	filter_set = list(set([f.filter for f in data]))
	print(filter_set)
	for filt in filter_set:
		filt_data = [x for x in data if x.filter == filt]
		min_mjd = min([x.mjd for x in filt_data])
		max_mjd = max([x.mjd for x in filt_data])
		print('Filter {}, min {} max {}'.format(filt, round(min_mjd, 2), round(max_mjd, 2)))
		while min_mjd < (max_mjd+day_thresh):
			
			mjd_set = [x for x in filt_data if x.mjd > (min_mjd-day_thresh) and x.mjd < (min_mjd+day_thresh)]
			print('len {}, range: {}-{}'.format(len(mjd_set), round(min_mjd-day_thresh, 2), round(min_mjd+day_thresh)))
			if len(mjd_set) > 1:
				mean_mag = np.mean([x.mag for x in mjd_set])
				mean_mag_err = np.mean([x.mag_err for x in mjd_set])
				std_mag = np.std([x.mag for x in mjd_set])
				clipped_mag = []

				for x in [f for f in mjd_set]:
					if abs(mean_mag - x.mag) < 1.5*std_mag:
						clipped_mag.append([x.mag, x.mag_err])

				if len(clipped_mag):
					mean_mag = np.mean([a[0] for a in clipped_mag])
					std_mag = np.std([a[0] for a in clipped_mag])
					mean_mag_err = np.mean([a[1] for a in clipped_mag])

					errprop = np.sqrt(std_mag**2 + mean_mag_err**2)
					mean_mjd = np.mean([x.mjd for x in mjd_set])

					mmd.append(
						exposureInfo(
							filt,
							mean_mag,
							errprop,
							mean_mjd,
							None,
							'charlie_data.csv'
						)
					)
			elif len(mjd_set) == 1:
				mmd.append(mjd_set[0])
			
			min_mjd = min_mjd + (day_thresh*2)

	return mmd


def mean_night_data(data):
	mnd = []
	filter_set = list(set([f.filter for f in data]))
	for i,filt in enumerate(filter_set):
		filt_data = [x for x in data if x.filter == filt]
		night_set = list(set([f.night for f in filt_data if f.night is not None]))
		for n in night_set:
			night_data = [f for f in filt_data if f.night == n]
			mean_mag = np.mean([x.mag for x in night_data])
			mean_mag_err = np.mean([x.mag_err for x in night_data])
			std_mag = np.std([x.mag for x in night_data])
			clipped_mag = []

			for x in [f for f in night_data]:
				if abs(mean_mag - x.mag) < 1.5*std_mag:
					clipped_mag.append([x.mag, x.mag_err])
			
			if len(clipped_mag):
				mean_mag = np.mean([a[0] for a in clipped_mag])
				std_mag = np.std([a[0] for a in clipped_mag])
				mean_mag_err = np.mean([a[1] for a in clipped_mag])


				errprop = np.sqrt(std_mag**2 + mean_mag_err**2)
				mean_mjd = np.mean([x.mjd for x in night_data])

				mnd.append(
					exposureInfo(
						filt,
						mean_mag,
						errprop,
						mean_mjd,
						n,
						'SN2021fxy_1m_phot_only.csv'
					)
				)

	return mnd


def find_lcogt_maxlightish():
	samfiles = lcfile('21fxy_phot.csv', 1, 5, 2, 3, 0, None, True)
	samdata = mean_night_data(samfiles.ingest())

	filters = ['g', 'r','i']
	for f in filters:
		filtdata = [x for x in samdata if x.filter == f]
		min_mag = min([x.mag for x in filtdata])
		min_utdate = [x.night for x in filtdata if x.mag == min_mag][0]
		print("maxlight {} utdate for filter {}: {}".format(min_mag, f, min_utdate))


def discrepancy():

	samfiles = lcfile('SN2021fxy_1m_phot_only.csv', 1, 5, 2, 3, 0, None, True)
	charliefiles = [
		lcfile('lcogt.dat', 0, 1, 2, 3),
		lcfile('lulin.dat', 0, 1, 2, 3),
		lcfile('thacher.dat', 0, 1, 2, 3)
	]

	#samdata = mean_night_data(samfiles.ingest())
	samdata = samfiles.ingest()
	charliedata = []
	for c in charliefiles:
		charliedata.extend(c.ingest())

	#charliedata = mean_mjd_data(charliedata)
	

	dataset = [['samdata', samdata], ['charliedata', charliedata]]

	universal_offset = 0
	xmin = 0
	offset_dict = {"u":-4.0, "B":-2.0, "g":-1.0, "V":0.0, "r":1.0, "i":1.5, "z":3.0}
	filters = ['g', 'r', 'i', 'z']
	#filters = ['z']
	colors = iter(cm.rainbow(np.linspace(0, 1, len(filters)*4)))

	for i,filt in enumerate(filters):
		print(filt)
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		font = {'family': 'serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14}

		ax1.set_xlabel(r'Days since $B_{max}$', fontdict=font)
		ax1.set_ylabel(r'Apparent Magnitude', fontdict=font)
		for f in dataset:
			data = f[1]
			filename = f[0]

			filt_data = [x for x in data if x.filter == filt]
			mag_dat = [x.mag for x in filt_data if x.mag < 21.5]
			#abmag_dat = [x.mag - d_mod for x in filt_data]
			mjd_dat = [x.mjd for x in filt_data if x.mag < 21.5]

			mjd_datNew = [x-universal_offset for x in mjd_dat]

			offset = offset_dict[filt]
			mag_datNew = [x + offset for x in mag_dat]

			if len(mjd_datNew):
				xmin = min(mjd_datNew) if min(mjd_datNew) < xmin else xmin

				strfset = str(offset) if offset < 0.0 else "+"+str(offset)
				sizes = [10 for x in mag_datNew]
				markerstyles = ['o' for x in mag_datNew]
				color = next(colors)
				
				#if filt == 'V':
				#	sizes[0] = 50
				#	markerstyles[0] = '*'

				for it in range(0, len(mjd_datNew), 1):
					ax1.scatter([mjd_datNew[it]], [mag_datNew[it]], s=[sizes[it]], color=color, marker=markerstyles[it])
				
				ax1.scatter([],[], label=filt+' ' +filename, color=color, marker=markerstyles[1], s=8)
		
		ax1.legend(fontsize=9)
		ax1.set_ylim(ax1.get_ylim()[::-1])
		print('plotting?')
		plt.savefig("discrep_{}.png".format(filt), format='png', dpi=1000)
		plt.show()

		
def main():
	lcfiles = [
		lcfile('SN2021fxy_1m_only_final.csv', 1, 5, 2, 3, 0, None, True),
		lcfile('dlt40_lc.txt', 0, -1, 2, 3, 1, None, True),
		lcfile('lcogt.dat', 0, 1, 2, 3),
		lcfile('lulin.dat', 0, 1, 2, 3),
		lcfile('thacher.dat', 0, 1, 2, 3)
	]

	overall_data = []
	for f in lcfiles:
		print(f.fname)
		if 'SN2021fxy_1m_phot_only.csv' in f.fname:
			#overall_data.extend(mean_night_data(f.ingest()))
			overall_data.extend(f.ingest())
		else:
			overall_data.extend(f.ingest())

	offset_dict = {"U":-5.0, "B":-4.0, "V":-3, 'clear':-2.2, "u":-2.0,  "g":-0.3, "r":0.2, "i":0.8, "z":1.6}
	filters = ["U", 'B', 'V', 'clear', 'u','g', 'r', 'i', 'z']

	dist = 33.88 #taken from :http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=NGC+5018
	d_mod = 5*math.log10(dist*(math.pow(10, 6))/10)
	print(d_mod)
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()

	xmin = 0
	universal_offset = 59305.130
	opticalspectradates = ['59347.4', '59327.4', '59313.6', '59308.1', '59306.1', '59304.9', '59300.7', '59297.9', '59297.5', '59295.5', '59292.6', '59291.3', '59291.1']

	ax1.plot((0, 0), (17, 8.3), "--k", alpha=0.5)

	for d in opticalspectradates:
		d = float(d) - universal_offset# + 2400000.5
		ax1.plot((d, d), (17.8,18.8), "-b")

	colors = iter(cm.rainbow(np.linspace(0, 1, len(filters))))

	for i,filt in enumerate(filters):
		filt_data = [x for x in overall_data if x.filter == filt]
		mag_dat = [x.mag for x in filt_data if x.mag < 21.5]
		abmag_dat = [x.mag - d_mod for x in filt_data]
		mjd_dat = [x.mjd for x in filt_data if x.mag < 21.5]

		mjd_datNew = [x-universal_offset for x in mjd_dat]

		offset = offset_dict[filt]
		mag_datNew = [x + offset for x in mag_dat]

		if len(mjd_datNew):
			xmin = min(mjd_datNew) if min(mjd_datNew) < xmin else xmin

			strfset = str(offset) if offset < 0.0 else "+"+str(offset)
			sizes = [10 for x in mag_datNew]
			markerstyles = ['o' for x in mag_datNew]
			color = next(colors)
			
			if filt == 'u':
				sizes = [13 for x in mag_datNew]
				#markerstyles[0] = '*'

			for it in range(0, len(mjd_datNew), 1):
				ax1.scatter([mjd_datNew[it]], [mag_datNew[it]], s=[sizes[it]], color=color, marker=markerstyles[it])
			
			ax1.scatter([],[], label=filt+strfset, color=color, marker=markerstyles[1], s=8)


	font = {'family': 'serif',
			'color':  'black',
			'weight': 'normal',
			'size': 14}

	ax1.set_xlabel(r'Days since $B_{max}$', fontdict=font)
	ax1.set_ylabel(r'Apparent Magnitude', fontdict=font)
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
	ax1.set_xlim(xmin-5, 65)
	print(xmin, 'earrrly')
	plt.savefig("SN2021fxy_lightcurve.png", format='png', dpi=1000)
	plt.show()

def write_snoopy():
	lcfiles = [
		lcfile('SN2021fxy_1m_only_final.csv', 1, 5, 2, 3, 0, None, True),
		lcfile('lcogt.dat', 0, 1, 2, 3),
		lcfile('lulin.dat', 0, 1, 2, 3),
		lcfile('thacher.dat', 0, 1, 2, 3)
	]

	overall_data = []
	for f in lcfiles:
		print(f.fname)
		if 'SN2021fxy_1m_only_final.csv' in f.fname:
			overall_data.extend(f.ingest())
		else:
			pass
			#overall_data.extend(f.ingest())
	output_for_snpy('SN2021fxy', 0.0094, 198.256560801, -19.5125038193, overall_data)


main()
#discrepancy()
#find_lcogt_maxlightish()
#write_snoopy()
#recarr = np.array([(filters), (mags), (mag_errs), (mjds), (nights)], dtype=dtype)
