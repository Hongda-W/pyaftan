import pyaftan
import matplotlib.pyplot as plt
import sys
import exceptions
import obspy
import timeit
"""
    OCT 2016, Hongda
"""
if len(sys.argv) != 4:
	print '!!! Usage: [Sac_file_name] [1sta_name] [2sta_name] !!!'
	raise ValueError
tr=obspy.read(sys.argv[1])[0]
atr=pyaftan.aftantrace(tr.data, tr.stats)
# atr1=atr.getpos() # atr1 is for the positive lag
# atr2=atr.getneg() # atr2 is for the negative lag
atr.makesym() # atr3 is for the mean of positive and negative lags
#start_time=timeit.default_timer()
# atr1.aftan(tmin=5., tmax=20., vmin=0.5, vmax=4.5, phvelname='/projects/howa1663/Code/pyaftan/ak135.disp') # apply aftan to all 3 traes
# atr2.aftan(tmin=5., tmax=20., vmin=0.5, vmax=4.5, phvelname='/projects/howa1663/Code/pyaftan/ak135.disp')
atr.aftan(tmin=5., tmax=40., vmin=0.5, vmax=4.5, phvelname='/projects/howa1663/Code/pyaftan/ak135.disp')
# atr1.get_snr()
# atr2.get_snr()
atr.get_snr()
# atr1.ftanparam.writeDISP(sys.argv[2]+'_'+sys.argv[3]+'_POS')
# atr2.ftanparam.writeDISP(sys.argv[2]+'_'+sys.argv[3]+'_NEG')
# atr.ftanparam.writeDISP(sys.argv[2]+'_'+sys.argv[3]+'_SYM')
#atr1.plotftan(plotflag=3, sacname=sys.argv[2]+'_'+sys.argv[3])
#print("Plotting first dispersion curve took %s seconds" % (timeit.default_timer()-start_time))
#plt.savefig(sys.argv[2]+'_'+sys.argv[3]+'_POS.png')
#plt.clf()
#atr2.plotftan(plotflag=3, sacname=sys.argv[2]+'_'+sys.argv[3])
#plt.savefig(sys.argv[2]+'_'+sys.argv[3]+'_NEG.png')
#plt.clf()
atr.plotftan(plotflag=3, sacname=sys.argv[2]+'_'+sys.argv[3])
plt.savefig(sys.argv[2]+'_'+sys.argv[3]+'_SYM.png')
#plt.clf()
