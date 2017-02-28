"""
    Plot the saved group dispersion curves, click to draw the reference dispersion curve, save the model
    left click to add data point, right click to delete data point
                    -- Hongda Wang; Feb 21, 2017 
"""
from matplotlib import pyplot as plt
from matplotlib import backend_bases as base
import numpy as np
import sys
import warnings
import exceptions
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
if len(sys.argv) != 2:
    print("!!!Usage: [Array name for data points]!!!")
    raise ValueError
scts = np.load(sys.argv[1]) # the scattered data points 
grp = sys.argv[1].split('.')[0]
MdlNam = "Md_" + grp # name for saving dispersion curve model

def PolyFit(data_x, data_y, fig, N=5):
    """
    Function using least squares to find ploynomial fit for the data. N is the degree of ploynomial that we use.
    """
    if len(data_x) != len(data_y):
        raise ValueError('Error: polynomial fit data must be of the same length')
    m = len(data_x)
    A = np.ones([m,N+1])
    for i in 1 + np.arange(N):
        A[:,i] = A[:,i-1] * data_x
    B = A.T.dot(A)
    y = A.T.dot(data_y)
    xs = np.linalg.solve(B,y)
    poly_fit = np.zeros(m)
    for i in np.arange(N+1):
        poly_fit += xs[i] * (data_x ** i)
    fig.plot(data_x, poly_fit, 'b-')

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.cid = line.figure.canvas.mpl_connect('key_press_event', self)
    def __call__(self, event):
        if isinstance(event, base.KeyEvent):
            if event.key == "enter": # Press 'enter' for polynomial fit
                FitCurve = PolyFit(self.xs, self.ys, N=7, fig=ax)
                self.line.figure.canvas.draw()
            if event.key == " ": # Press 'space' for savgol_filter smooth
                xs = np.asarray(self.xs)
                ys = np.asarray(self.ys)
                itp = interp1d(xs, ys, kind='linear')
                xx = np.linspace(xs.min(), xs.max(), 1000)
                ys_sg = savgol_filter(itp(xx), 101, 3)
                ax.plot(xx, ys_sg, 'y-')
                self.line.figure.canvas.draw()
        if isinstance(event, base.MouseEvent):
            print('click', event)
            if event.inaxes!=self.line.axes: return
            # print("You pressed: "+str(event.button))
            if event.button == 3: # right click to remove last point
                self.xs = self.xs[:-1]
                self.ys = self.ys[:-1]
            if event.button == 1: # left click to add point
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
            self.line.set_data(self.xs, self.ys)
            self.line.set_color("red")
            self.line.figure.canvas.draw()
print("Left click to add data point, right click to delete last data point, close window to finish!")
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Click to draw reference dispersion curve for '+ grp)
ax.plot(scts[0,:], scts[1,:], 'c.', markersize=0.4)
ax.set_xlabel("Period (s)")
ax.set_ylabel("Velocity (km/s)")
ax.set_xlim([0, 45])
ax.set_ylim([0, 5])

line, = ax.plot([],[])  # empty line
linebuilder = LineBuilder(line)
plt.show()
model=np.vstack([linebuilder.xs,linebuilder.ys])
np.save(MdlNam, model)
