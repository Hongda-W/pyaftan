import FTAN_4_JdF
import numpy as np
import matplotlib.pyplot as plt
import sys
import exceptions
import obspy
import glob
import re
"""
    run FTAN_4_JdF for all station pairs, save their dispersion result
    Feb 2017, Hongda
"""
if len(sys.argv) != 2:
	print '!!! Usage: [Directory for stacked cross-correlation] !!!'
	raise ValueError
with open('/projects/howa1663/Code/NoisePy/station_dep_4_plot.lst') as f:
    lines = f.readlines()
# [m.group(0) for l in lines for m in [regex.search(l)] if m]
cc = np.array([[],[]])
cs = np.array([[],[]])
cd = np.array([[],[]])
ss = np.array([[],[]])
sd = np.array([[],[]])
dd = np.array([[],[]])
for SAC_file in glob.glob(sys.argv[1]+'/*/*.SAC'):
    staNam1 = SAC_file.split('/')[-1].split('_')[1] # station name, specified for '*/COR_FS09D_M17D.SAC' files
    staNam2 = SAC_file.split('/')[-1].split('_')[-1].split('.')[0]
    pat1 = re.compile(staNam1)
    match1 = filter(pat1.match, lines)
    pat2 = re.compile(staNam2)
    match2 = filter(pat2.match, lines)
    if len(match1) <= 2: # ignore multiple stations with the same staName
        flag1 = int(match1[0].split()[-1])
        if flag1 == 3: # change the flag for deep water stations, in order to get unique flag for different pairs
            flag1 = 4
    if len(match2) <= 2:
        flag2 = int(match2[0].split()[-1])
        if flag2 == 3:
            flag2 = 4
    print("Start FTAN for "+staNam1+"_"+staNam2)
    flag = flag1 + flag2
    tr=obspy.core.read(SAC_file)[0]
    atr=FTAN_4_JdF.aftantrace(tr.data, tr.stats)
    atr.makesym()
    atr.aftan(tmin=1., tmax=40., vmin=0.2, vmax=4.5)
    if atr.ftanparam.nfout1_1 == 0:
        print("No ftan result for "+staNam1+"_"+staNam2)
        continue
    atr.get_snr()
    per = atr.ftanparam.arr1_1[1,:atr.ftanparam.nfout1_1]
    gvel = atr.ftanparam.arr1_1[2,:atr.ftanparam.nfout1_1]
    TV_pair = np.vstack([per,gvel])
    if flag == 2: # cont-cont
        cc = np.append(cc,TV_pair,axis=1)
    elif flag == 3: # cont-shallow
        cs = np.append(cs,TV_pair,axis=1)
    elif flag == 5: # cont-deep
        cd = np.append(cd,TV_pair,axis=1)
    elif flag == 4: # shallow-shallow
        ss = np.append(ss,TV_pair,axis=1)
    elif flag == 6: # shallow-deep
        sd = np.append(sd,TV_pair,axis=1)
    elif flag == 8: # deep-deep
        dd = np.append(dd,TV_pair,axis=1)
    else:
        print( "Wrong value for flag of station water depth type!!! flag1 is "+str(flag1)+" flag2 is "+str(flag2) )
        raise ValueError
    #atr.ftanparam.writeDISP('/lustre/janus_scratch/howa1663/CC_JdF/DISP_4_GrpV/'+staNam1+'_'+staNam2+'_SYM')
np.save('ForAll_cc', cc)
np.save('ForAll_cs', cs)
np.save('ForAll_cd', cd)
np.save('ForAll_ss', ss)
np.save('ForAll_sd', sd)
np.save('ForAll_dd', dd)
