import pyaftan
import obspy
import matplotlib.pyplot as plt
import numpy as np
import timeit



# tr=obspy.read('/lustre/janus_scratch/howa1663/CC_JdF/Stack/J27B/COR_J27B_J48B.SAC')[0]
tr=obspy.read('/lustre/janus_scratch/howa1663/CC_JdF/Stack/J33C/COR_J33C_M08C.SAC')[0]
atr1=pyaftan.aftantrace(tr.data, tr.stats)
# atr2=pyaftan.aftantrace(tr.data, tr.stats)
# aftan analysis using pyaftan
# t1=timeit.default_timer()
atr1.makesym()
atr1.aftan(tmin=5., tmax=50., vmin=0.5, vmax=5., phvelname='ak135.disp')
atr1.get_snr()
atr1.ftanparam.writeDISP('J27B_J48B')
# for i in xrange(10):
#     print i
#     atr1.aftan(tmin=5., tmax=20., vmin=2.5, vmax=3.5, phvelname='ak135.disp')
# c = f(10000, 10000)
# t2=timeit.default_timer()
# print t2-t1
atr1.plotftan(plotflag=3)
# # plt.suptitle('pyaftan results')
# # aftan analysis using compiled fortran library
# # atr2.aftanf77(tmin=5., tmax=30., vmin=2.5, vmax=4.5, phvelname='ak135.disp')
# # atr2.plotftan(plotflag=3)
# # plt.suptitle('fortran77 aftan results')
plt.show()
