import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
try:
    plt.style.use('ggplot')
except:
    print("ggplot Not available")
def svgl_smth(x, y):
    itp = interp1d(x, y, kind='linear')
    xx = np.linspace(x.min(), x.max(), 1000)
    ys_sg = savgol_filter(itp(xx), 151, 3)
    return ys_sg, xx
cc = np.load("Md_ForAll_cc.npy")
cs = np.load("Md_ForAll_cs.npy")
cd = np.load("Md_ForAll_cd.npy")
ss = np.load("Md_ForAll_ss.npy")
sd = np.load("Md_ForAll_sd.npy")
dd = np.load("Md_ForAll_dd.npy")
cc_sg, cc_xx = svgl_smth(cc[0,:], cc[1,:])
cs_sg, cs_xx = svgl_smth(cs[0,:], cs[1,:])
cd_sg, cd_xx = svgl_smth(cd[0,:], cd[1,:])
ss_sg, ss_xx = svgl_smth(ss[0,:], ss[1,:])
sd_sg, sd_xx = svgl_smth(sd[0,:], sd[1,:])
dd_sg, dd_xx = svgl_smth(dd[0,:], dd[1,:])

plt.figure()
plt.plot(cc[0,:], cc[1,:], color='r', label='cc', linewidth=2)
plt.plot(cs[0,:], cs[1,:], color='y', label='cs', linewidth=2)
plt.plot(cd[0,:], cd[1,:], color='b', label='cd', linewidth=2)
plt.plot(ss[0,:], ss[1,:], color='c', label='ss', linewidth=2)
plt.plot(sd[0,:], sd[1,:], color='m', label='sd', linewidth=2)
plt.plot(dd[0,:], dd[1,:], color='g', label='dd', linewidth=2)
plt.legend(loc=4)
plt.xlim(0,40)
plt.xlabel("Period (sec)")
plt.ylabel("Group Velocity (km/sec)")

plt.figure()
plt.plot(cc_xx, cc_sg, color='r', label='cc', linewidth=2)
plt.plot(cs_xx, cs_sg, color='y', label='cs', linewidth=2)
plt.plot(cd_xx, cd_sg, color='b', label='cd', linewidth=2)
plt.plot(ss_xx, ss_sg, color='c', label='ss', linewidth=2)
plt.plot(sd_xx, sd_sg, color='m', label='sd', linewidth=2)
plt.plot(dd_xx, dd_sg, color='g', label='dd', linewidth=2)
plt.legend(loc=4)
plt.xlim(0,40)
plt.xlabel("Period (sec)")
plt.ylabel("Group Velocity (km/sec)")
plt.show()
