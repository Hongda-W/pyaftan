import matplotlib.pyplot as plt
import numpy as np
import sys
if len(sys.argv) != 2:
    print("!!!Usage: [Array name for data points]!!!")
scts = np.load(sys.argv[1]) # load data points
# scts = scts[scts[:,0].argsort()] # sort by periods
bins = np.linspace(1, 42, 10001)
boxes = np.arange(0, 5, 0.01)
h = 0.041
hh = 0.01
per = scts[0, :]
vel = scts[1, :]
cnts = np.array([])
mdlV = np.array([])
for i in bins:
    inds = np.logical_and(per > (i-h/2), per <= (i+h/2))
    vs = vel[inds]
    cnts = np.array([])
    for j in boxes:
        iind = np.logical_and(vs > j,  vs <= j+hh ) * 1
        cnts = np.append(cnts, iind.sum())
    lclM = np.where(cnts==cnts.max())[0]
    if len(lclM) > 1:
        print("Find max counts in multiple boxes")
    usefulV = vs[np.logical_and(vs > boxes[lclM[-1]],  vs <= boxes[lclM[-1]]+hh )]
    mdlV = np.append(mdlV, usefulV.mean())
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(scts[0,:], scts[1,:], 'c.', markersize=0.4)
ax.set_xlabel("Period (s)")
ax.set_ylabel("Velocity (km/s)")
ax.set_xlim([0, 45])
ax.set_ylim([0, 5])
ax.plot(bins, mdlV, 'r-')
plt.show()
