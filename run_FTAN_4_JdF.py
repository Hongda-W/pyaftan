import FTAN_4_JdF
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
atr=FTAN_4_JdF.aftantrace(tr.data, tr.stats)
# atr1=atr.getpos() # atr1 is for the positive lag
# atr2=atr.getneg() # atr2 is for the negative lag
atr.makesym()
atr.aftan(tmin=1., tmax=40., vmin=0.2, vmax=4.5)
atr.get_snr()
#atr.ftanparam.writeDISP('/lustre/janus_scratch/howa1663/CC_JdF/DISP_4_GrpV'+sys.argv[2]+'_'+sys.argv[3]+'_SYM')
atr.plot_group(title=sys.argv[2]+'_'+sys.argv[3])
plt.show()
# plt.savefig('Group_'+sys.argv[2]+'_'+sys.argv[3]+'_SYM.png')
#plt.clf()
