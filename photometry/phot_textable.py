import numpy as np
import matplotlib.pyplot as plt
import math

class exposureInfo:
	def __init__(self, f, m, m_e, mjd, n):
		self.filter = f
		self.mag = m
		self.mag_err = m_e
		self.mjd = mjd
		self.night = n

class nightobj:
	def __init__(self, date, mjd=None):
		self.date = date
		self.mjd = mjd
		self.filtdict = {}

class filterinfo:
	def __init__(self, mag, mag_err):
		self.mag = mag
		self.mag_err = mag_err

def maketable(filters, outfile):
	filt_columns = ""
	filt_columns_names = ""
	filt_columns_orient = ""
	for f in filters:
		filt_columns_orient += "r"
		filt_columns += "c"
		filt_columns_names += " & \\colhead{"+f+"}"

	#\begin{deluxetable*}{lrr}%[!t]
	#\tablecaption{Properties of Col~I and Tri~II \label{tab:params}}
	#\tablehead{\colhead{Parameter} & \colhead{Col~I} & \colhead{Tri~II}}
	#\startdata

	out = "\\begin{deluxetable*}{ccc"+filt_columns+"}\n"
	out += "\\tablecaption{Photometric data of CSP15AAK \\label{tab:params}}\n"
	out += "\\tablehead{\\colhead{Date} & \\colhead{JD} & \\colhead{Phase}"+filt_columns_names+"}\n"
	out += "\\startdata\n"

	for n in nights:
		if any(x in n.filtdict.keys() for x in filters):
			out += "{} & {} & {}".format(n.date, n.mjd, "")
			for f in filters:
				if f in n.filtdict.keys():
					out += " & ${} \\pm {}$".format(n.filtdict[f].mag, n.filtdict[f].mag_err)
				else:
					out == " &"
			out += " \\\\\n"

	out += "\\enddata\n"
	out += "\\end{deluxetable*}"

	of = open(outfile, "w")
	of.write(out)
	of.close()

lcFile = "CSP15aak_phot_csp.dat"
filters = ['u','B','g','V','r','i','Y','J','H']
filters_opt = ['u','B','g','V','r','i']
filters_nir = ['Y','J','H'] 
lcFile_o = open(lcFile).readlines()[1:]

data = []
for row in lcFile_o:
	rs = row.split()
	data.append(exposureInfo(rs[3], float(rs[5]),float(rs[6]),float(rs[9]),rs[0]))


dates = sorted(list(set([x.night for x in data])))

nights = []
for date in dates:
	night_inst = nightobj(date)
	for d in data:
		if d.night == date:
			night_inst.mjd = d.mjd
			if d.filter not in night_inst.filtdict.keys():
				night_inst.filtdict[d.filter] = filterinfo(d.mag, d.mag_err)
			else:
				print "Already in nightly exposure: ", d.night, d.filter
	nights.append(night_inst)
print len(nights)


maketable(filters_opt, "phot_tex_table_opt.txt")
maketable(filters_nir, "phot_tex_table_nir.txt")



