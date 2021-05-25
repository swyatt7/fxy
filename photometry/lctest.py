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

lcFile = "CSP15aak_phot_csp.dat"
lcFile_o = open(lcFile).readlines()[1:]

data = []
for row in lcFile_o:
	rs = row.split()
	data.append(exposureInfo(rs[3], float(rs[5]),float(rs[6]),float(rs[9]),rs[0]))

bdat = [x for x in data if x.filter == 'B']
bmags = [x.mag for x in bdat]
mjds = [x.mjd for x in bdat]

peakB_i = bmags.index(min(bmags))
offset = mjds[peakB_i]
mjds = [int(x - offset) for x in mjds]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.scatter(mjds,bmags)
ax1.set_xlabel("Days since max B magnitude")
ax1.set_ylabel("Apparent Magnitude")

ax2.set_ylabel("Absolute Magnitude")


ax1Ys = ax1.get_yticks()
ax2Ys = []
dist = 22.650
d_mod = 5*math.log10(dist*(math.pow(10, 6))/10)

for X in ax1Ys:
    ax2Ys.append(int(X - d_mod))

ax2.set_yticks(ax1Ys)
ax2.set_ybound(ax1.get_ybound())
ax2.set_yticklabels(ax2Ys)
ax1.set_ylim(ax1.get_ylim()[::-1])
ax2.set_ylim(ax2.get_ylim()[::-1])
#title = ax1.set_title("Upper x-axis ticks are lower x-axis ticks doubled!")
#title.set_y(1.1)
#fig.subplots_adjust(top=0.85)
plt.show()