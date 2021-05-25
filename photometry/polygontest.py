import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


fig, ax = plt.subplots()
patches = []

nuvredx = [-15,-14,-12,-8,0,12,19,-5,-10,-15]
nuvredy = [1.2,1.3,1.1,0.1,0.3,1.4,1.4,-0.6,-0.2,1.2]
nuvred = []
for i,x in enumerate(nuvredx):
	nuvred.append([x, nuvredy[i]])

nuvbluex = [-16,-16,-13,-10,-6,0,20,20,-3,-8,-11,-16]
nuvbluey = [0.5,0.7,0.5,-0.5,-0.73,-0.3,1.2,0.7,-0.95,-0.75,-0.6,0.5]
nuvblue = []

for i,x in enumerate(nuvbluex):
	nuvblue.append([x, nuvbluey[i]])

patches.append(Polygon(nuvred, color='red'))
patches.append(Polygon(nuvblue, color='blue'))

colors = [55, 15]
print colors
p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
p.set_array(np.array(colors))
p.set_clim([15, 55])
ax.add_collection(p, autolim=True)
ax.autoscale_view()
#fig.colorbar(p, ax=ax)

plt.show()
