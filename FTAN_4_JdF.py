# -*- coding: utf-8 -*-
"""
A python module for applying aftan analysis to the ambient noise cross-correlation data

A revised verison from Lili's pyaftan.py

                                                                   -- Hongda Feb, 2017
"""

import obspy
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import os
import warnings
from scipy.signal import argrelmax, argrelmin, argrelextrema
import scipy.interpolate
try:
    import pyfftw
    useFFTW=True
except:
    useFFTW=False
try:
    import aftan
    isaftanf77=True
except:
    isaftanf77=False


def _aftan_gaussian_filter(alpha, omega0, ns, indata, omsArr):
    """Internal Gaussian filter function used for aftan
    """
    om2 = -(omsArr-omega0)*(omsArr-omega0)*alpha/omega0/omega0
    b=np.exp(om2)
    b[np.abs(om2)>=40.]=0.
    filterred_data=indata*b
    filterred_data[ns/2:]=0
    filterred_data[0]/=2
    filterred_data[ns/2-1]=filterred_data[ns/2-1].real+0.j
    return filterred_data

class ftanParam(object):
    """ An object to handle ftan output parameters
    ===========================================================================
    Basic FTAN parameters:
    nfout1_1 - output number of frequencies for arr1, (integer*4)
    arr1_1   - preliminary results.
                Description: real*8 arr1(8,n), n >= nfin)
                arr1_1[0,:] -  central periods, s
                arr1_1[1,:] -  observed periods, s
                arr1_1[2,:] -  group velocities, km/s
                arr1_1[3,:] -  phase velocities, km/s or phase if nphpr=0, rad
                arr1_1[4,:] -  amplitudes, Db
                arr1_1[5,:] -  discrimination function
                arr1_1[6,:] -  signal/noise ratio, Db
                arr1_1[7,:] -  maximum half width, s
                arr1_1[8,:] -  amplitudes
    ===========================================================================
    """
    def __init__(self):
        # Parameters for first iteration
        self.nfout1_1=0
        self.arr1_1=np.array([])
        self.tamp_1=0.
        self.nrow_1=0
        self.ncol_1=0
        self.ampo_1=np.array([],dtype='float32')
        self.ierr_1=0
        # Parameters for second iteration
        self.nfout1_2=0
        self.arr1_2=np.array([])
        self.nfout2_2=0
        self.arr2_2=np.array([])
        self.tamp_2=0.
        self.nrow_2=0
        self.ncol_2=0
        self.ampo_2=np.array([])
        self.ierr_2=0
        # Flag for existence of predicted phase dispersion curve
        self.preflag=False
        self.station_id=None

    def writeDISP(self, fnamePR):
        """
        Write FTAN parameters to DISP files given a prefix.
        fnamePR: file name prefix
        _1_DISP.0: arr1_1
        _1_DISP.1: arr2_1
        _2_DISP.0: arr1_2
        _2_DISP.1: arr2_2
        """
        if self.nfout1_1!=0:
            f10=fnamePR+'_1_DISP.0'
            Lf10=self.nfout1_1
            outArrf10=np.arange(Lf10)
            for i in xrange(9):
                outArrf10=np.append(outArrf10, self.arr1_1[i,:Lf10])
            outArrf10=outArrf10.reshape((10,Lf10))
            outArrf10=outArrf10.T
            np.savetxt(f10, outArrf10, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf %12.4lf %12.4lf')
        return
    
    def writeDISPbinary(self, fnamePR):
        """
        Write FTAN parameters to DISP files given a prefix.
        fnamePR: file name prefix
        _1_DISP.0: arr1_1
        _1_DISP.1: arr2_1
        _2_DISP.0: arr1_2
        _2_DISP.1: arr2_2
        """
        f10=fnamePR+'_1_DISP.0'
        np.savez(f10, self.arr1_1, np.array([self.nfout1_1]) )
        f11=fnamePR+'_1_DISP.1'
        np.savez(f11, self.arr2_1, np.array([self.nfout2_1]) )
        f20=fnamePR+'_2_DISP.0'
        np.savez(f20, self.arr1_2, np.array([self.nfout1_2]) )
        f21=fnamePR+'_2_DISP.1'
        np.savez(f21, self.arr2_2, np.array([self.nfout2_2]) )
        return
    

    def FTANcomp(self, inftanparam, compflag=1):
        """
        Compare aftan results for two ftanParam objects.
        """
        fparam1=self
        fparam2=inftanparam
        if compflag==1:
            obper1=fparam1.arr1_1[1,:fparam1.nfout1_1]
            gvel1=fparam1.arr1_1[2,:fparam1.nfout1_1]
            phvel1=fparam1.arr1_1[3,:fparam1.nfout1_1]
            obper2=fparam2.arr1_1[1,:fparam2.nfout1_1]
            gvel2=fparam2.arr1_1[2,:fparam2.nfout1_1]
            phvel2=fparam2.arr1_1[3,:fparam2.nfout1_1]
        elif compflag==2:
            obper1=fparam1.arr2_1[1,:fparam1.nfout2_1]
            gvel1=fparam1.arr2_1[2,:fparam1.nfout2_1]
            phvel1=fparam1.arr2_1[3,:fparam1.nfout2_1]
            obper2=fparam2.arr2_1[1,:fparam2.nfout2_1]
            gvel2=fparam2.arr2_1[2,:fparam2.nfout2_1]
            phvel2=fparam2.arr2_1[3,:fparam2.nfout2_1]
        elif compflag==3:
            obper1=fparam1.arr1_2[1,:fparam1.nfout1_2]
            gvel1=fparam1.arr1_2[2,:fparam1.nfout1_2]
            phvel1=fparam1.arr1_2[3,:fparam1.nfout1_2]
            obper2=fparam2.arr1_2[1,:fparam2.nfout1_2]
            gvel2=fparam2.arr1_2[2,:fparam2.nfout1_2]
            phvel2=fparam2.arr1_2[3,:fparam2.nfout1_2]
        else:
            obper1=fparam1.arr2_2[1,:fparam1.nfout2_2]
            gvel1=fparam1.arr2_2[2,:fparam1.nfout2_2]
            phvel1=fparam1.arr2_2[3,:fparam1.nfout2_2]
            obper2=fparam2.arr2_2[1,:fparam2.nfout2_2]
            gvel2=fparam2.arr2_2[2,:fparam2.nfout2_2]
            phvel2=fparam2.arr2_2[3,:fparam2.nfout2_2]
        plb.figure()
        ax = plt.subplot()
        ax.plot(obper1, gvel1, '--k', lw=3) #
        ax.plot(obper2, gvel2, 'bo', markersize=5)
        plt.xlabel('Period(s)')
        plt.ylabel('Velocity(km/s)')
        plt.title('Group Velocity Comparison')
        if (fparam1.preflag and fparam2.preflag):
            plb.figure()
            ax = plt.subplot()
            ax.plot(obper1, phvel1, '--k', lw=3) #
            ax.plot(obper2, phvel2, 'bo', markersize=5)
            plt.xlabel('Period(s)')
            plt.ylabel('Velocity(km/s)')
            plt.title('Phase Velocity Comparison')
        return

class InputFtanParam(object): ###
    """
    A subclass to store input parameters for aftan analysis and SNR Analysis
    ===============================================================================================================
    Parameters:
    pmf         - flag for Phase-Matched-Filtered output (default: Fasle)
    piover4     - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
    vmin        - minimal group velocity, km/s
    vmax        - maximal group velocity, km/s
    tmin        - minimal period, s
    tmax        - maximal period, s
    tresh       - treshold for jump detection, usualy = 10, need modifications
    ffact       - factor to automatic filter parameter, usualy =1
    taperl      - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
    snr         - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
    fmatch      - factor to length of phase matching window
    fhlen       - half length of Gaussian width
    nfin        - number of initial period points
    npoints     - number of continuous points in jump correction
    perc        - output segment
    predV       - predicted phase velocity curve, period = predV[:, 0],  Vph = predV[:, 1]
    ===============================================================================================================
    """
    def __init__(self):
        self.pmf     = True
        self.piover4 = -1.0
        self.vmin    = 1.5
        self.vmax    = 5.0
        self.tmin    = 4.0
        self.tmax    = 70.0
        self.tresh   = 20.0
        self.ffact   = 1.0
        self.taperl  = 1.0
        self.snr     = 0.2
        self.fmatch  = 1.0
        self.fhlen   = 0.008
        self.nfin    = 64
        self.npoints = 3
        self.perc    = 50
        self.predV   = np.array([])
    
class aftantrace(obspy.core.trace.Trace):
    """
    aftantrace:
    A derived class inherited from obspy.core.trace.Trace. This derived class have a variety of new member functions
    """
    def init_ftanParam(self):
        """Initialize ftan parameters
        """
        self.ftanparam=ftanParam()
        
    def reverse(self):
        """Reverse the trace
        """
        self.data=self.data[::-1]
        return
    
    def makesym(self):
        """Turn the double lagged cross-correlation data to one single lag
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!')
        if self.stats.npts%2!=1:
            raise ValueError('Error: Incompatible begin and end time!')
        nhalf=(self.stats.npts-1)/2+1
        neg=self.data[:nhalf]
        pos=self.data[nhalf-1:self.stats.npts]
        neg=neg[::-1]
        self.data= (pos+neg)/2 
        self.stats.npts=nhalf
        self.stats.starttime=self.stats.starttime+self.stats.sac.e
        self.stats.sac.b=0.
        return

    def getneg(self):
        """Get the negative lag of a cross-correlation record
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!')
        negTr=self.copy()
        t=self.stats.starttime
        L=(int)((self.stats.npts-1)/2)+1
        negTr.data=negTr.data[:L]
        negTr.data=negTr.data[::-1]
        negTr.stats.npts=L
        negTr.stats.sac.b=0.
        negTr.stats.starttime=t-self.stats.sac.b
        return negTr

    def getpos(self):
        """Get the positive lag of a cross-correlation record
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!')
        posTr=self.copy()
        t=self.stats.starttime
        L=(int)((self.stats.npts-1)/2)+1
        posTr.data=posTr.data[L-1:]
        posTr.stats.npts=L
        posTr.stats.sac.b=0.
        posTr.stats.starttime=t-self.stats.sac.b
        return posTr
    
    def aftan(self, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, tmax=30.0, tresh=20.0, ffact=1.0,
            taperl=1.0, snr=0.2, fmatch=1.0, nfin=64, perc=50., phvelname='', predV=np.array([]), grpmdl=''):
        """ (Automatic Frequency-Time ANalysis) aftan analysis:
        ===========================================================================================================
        Input Parameters:
        piover4    - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin       - minimal group velocity, km/s
        vmax       - maximal group velocity, km/s
        tmin       - minimal period, s
        tmax       - maximal period, s
        tresh      - threshold for jump detection, usualy = 10, need modifications
        ffact      - factor to automatic filter parameter, usualy =1
        taperl     - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        snr        - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
        fmatch     - factor to length of phase matching window
        nfin       - number of initial period points
        perc       - output segment
        phvelname  - predicted phase velocity file name
        predV      - predicted phase velocity curve, period = predV[:, 0],  Vph = predV[:, 1]
        grpnam     - group velocity model file name (Added by Hongda Feb 2017)
        
        Output:
        self.ftanparam, a object of ftanParam class, to store output aftan results
        ===========================================================================================================
        """
        # preparing for data
        try:
            self.ftanparam
        except:
            self.init_ftanParam()
        try:
            dist=self.stats.sac.dist
        except:
            dist, az, baz=obspy.geodetics.base.gps2dist_azimuth(self.stats.sac.evla, self.stats.sac.evlo,
                                self.stats.sac.stla, self.stats.sac.stlo) # distance is in m
            self.stats.sac.dist=dist/1000.
            dist=dist/1000.
        if predV.size != 0:
            self.ftanparam.preflag=True
        elif os.path.isfile(phvelname):
            predV=np.loadtxt(phvelname)
            self.ftanparam.preflag=True
        else:
            warnings.warn('No predicted dispersion curve for:'+self.stats.network+'.'+self.stats.station, UserWarning, stacklevel=1)
        # Basic aftan
        grpVmdl = np.load(grpnam)
        if dist > 10**(-1):
            self._aftanpg(piover4=piover4, vmin=vmin, vmax=vmax, tmin=tmin, tmax=tmax, tresh=tresh, ffact=ffact, taperl=taperl,
                nfin=nfin, perc=perc, predV=predV, grpVmdl=grpVmdl)
        else:
            print("Interstation distance is too small, skip FTAN")
    def _aftanpg(self, piover4, vmin, vmax, tmin, tmax, tresh, ffact, taperl, nfin, perc, predV, grpVmdl):
        """ Basic aftan analysis, internal function
        ===========================================================================================================
        Input Parameters:
        piover4    - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin       - minimal group velocity, km/s
        vmax       - maximal group velocity, km/s
        tmin       - minimal period, s
        tmax       - maximal period, s
        tresh      - treshold for jump detection, usualy = 10, need modifications
        ffact      - factor to automatic filter parameter, usualy =1
        taperl     - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        nfin       - number of initial period points
        npoints    - number of continuous points in jump correction
        perc       - output segment
        predV      - predicted phase velocity curve, period = predV[:, 0],  Vph = predV[:, 1]
        grpVmdl    - group velocity model, grpVmdl[0,:]: periods; grpVmdl[1,:]: group velocity
        ===========================================================================================================
        Output:
        self.ftanparam, a object of ftanParam class, to store output aftan results
        """
        if self.ftanparam.preflag:
            phprper=predV[:,0]
            phprvel=predV[:,1]
            nprpv = predV[:,0].size
        else:
            nprpv=0
            phprper=np.array([])
            phprvel=np.array([])
        dt=self.stats.delta
        tb=self.stats.sac.b
        nsam=self.stats.npts
        dist=self.stats.sac.dist
        # alpha=ffact*20.*np.sqrt(dist/1000.)
        alpha=ffact*20.
        # number of samples for tapering, left and right end
        ntapb = int(round(taperl*tmax/dt))
        ntape = int(round(tmax/dt))
        omb = 2.0*np.pi/tmax
        ome = 2.0*np.pi/tmin
        # tapering seismogram
        nb = int(max(2, round((dist/vmax-tb)/dt)))
        tamp = (nb-1)*dt+tb
        ne = int(min(nsam, round((dist/vmin-tb)/dt)))
        nrow = nfin
        ncol = ne-nb+1
        tArr=np.arange(ne-nb+1)*dt+tb
        tArr[tArr==0.]=-1.
        vArr=dist/tArr
        tdata, ncorr=self.taper( max(nb, ntapb+1), min(ne, self.stats.npts-ntape), ntapb, ntape)
        # prepare for FFT
        ns=max(1<<(ncorr-1).bit_length(), 2**12)  # different !!!
        domega = 2.*np.pi/ns/dt
        step =(np.log(omb)-np.log(ome))/(nfin -1)
        omegaArr=np.exp(np.log(ome)+np.arange(nfin)*step)
        perArr=2.*np.pi/omegaArr
        # FFT
        if useFFTW:
            fftdata=pyfftw.interfaces.numpy_fft.fft(tdata, ns)
        else:
            fftdata=np.fft.fft(tdata, ns)
        omsArr=np.arange(ns)*domega
        phaArr=np.zeros((ne+3-nb, nfin))
        ampo=np.zeros((ne+3-nb, nfin))
        amp=np.zeros((ne+3-nb, nfin))
        #  main loop by frequency
        for k in xrange(nfin):
            # Gaussian filter
            filterS=_aftan_gaussian_filter(alpha=alpha, omega0=omegaArr[k], ns=ns, indata=fftdata, omsArr=omsArr)
            # return fs, filterS
            if useFFTW:
                filterT=pyfftw.interfaces.numpy_fft.ifft(filterS, ns)
            else:
                filterT=np.fft.ifft(filterS, ns)
            # need to multiply by 2 due to zero padding of negative frequencies
            # but NO NEED to divide by ns due to the difference of numpy style and FFTW style
            filterT=2.*filterT 
            phaArr[:,k]=np.arctan2(np.imag(filterT[nb-2:ne+1]), np.real(filterT[nb-2:ne+1]))
            ampo[:,k] = np.abs(filterT[nb-2:ne+1])
            amp[:,k] = 20.*np.log10(ampo[:,k])
        # normalization amp diagram to 100 Db with three decade cutting
        amax=amp.max()
        amp= amp+100.-amax
        amp[amp<40.]=40.
        tim1=np.zeros(nfin)
        tvis1=np.zeros(nfin)
        ampgr1=np.zeros(nfin)
        grvel1=np.zeros(nfin)
        snr1=np.zeros(nfin)
        wdth1=np.zeros(nfin)
        phgr1=np.zeros(nfin)
        ind_all=[]
        ipar_all=[]
        for k in xrange(nfin):
            ampk=amp[:,k]
            ampok=ampo[:,k]
            ind_localmax=argrelmax(ampk)[0] # index for local maxima
            omega=omegaArr[k]
            if ind_localmax.size==0:
                ind_localmax=np.append(ind_localmax, ne+1-nb)
            ind_all.append(ind_localmax)
            dph, tm, ph, t=self._fmax(amp=ampk, pha=phaArr[:,k], ind=ind_localmax, om=omega, piover4=piover4)
            imax=tm.argmax()
            ipar=np.zeros((6, ind_localmax.size))
            ipar[0, :]  = (nb+ind_localmax-2.+t)*dt # note the difference with aftanf77, due to ind_localmax
            ipar[1, :]  = 2*np.pi*dt/dph
            ipar[2, :]  = tm
            ipar[5, :]  = ph
            if ind_localmax.size==1:
                lmindex=(ampok[:ind_localmax[0]]).argmin()
                rmindex=(ampok[ind_localmax[0]:]).argmin()
                lm=(ampok[:ind_localmax[0]])[lmindex]
                rm=(ampok[ind_localmax[0]:])[rmindex]
            else:
                splitArr=np.split(ampok, ind_localmax)
                minArr=np.array([])
                minindexArr=np.array([])
                for tempArr in splitArr:
                    temp_ind_min=tempArr.argmin()
                    minArr=np.append(minArr, tempArr[temp_ind_min])
                    minindexArr=np.append(minindexArr, temp_ind_min)
                lm=minArr[:-1]
                rm=minArr[1:]
                minindexArr[1:]=minindexArr[1:]+ind_localmax
                lmindex=minindexArr[:-1]
                rmindex=minindexArr[1:]
            ipar[3,:] = 20.*np.log10(ampok[ind_localmax]/np.sqrt(lm*rm))
            ipar[4,:] = (np.abs(ind_localmax-lmindex)+np.abs(ind_localmax-rmindex))/2.*dt
            tim1[k]   = ipar[0,imax]
            tvis1[k]  = ipar[1,imax]
            ampgr1[k] = ipar[2,imax]
            if tvis1[k] <= (dist / 12): # delta/12 dispersion measurement limit
                grvel1[k] = dist/(tim1[k] +tb)
            else:
                grvel1[k] = 0
            snr1[k]   = ipar[3,imax]
            wdth1[k]  = ipar[4,imax] ### Note half width is not completely the same as ftanf77, need further check!!!
            phgr1[k]  = ipar[5,imax]
            ipar_all.append(ipar)
        nfout1=nfin
        per1=perArr
        shape=per1.shape # do not need ftrig1 set it to zero -hongda
        ftrig1=np.zeros(shape) #hongda
        if nprpv!=0:
            phV1=self._phtovel(per=tvis1, U=grvel1, pha=phgr1, npr=nprpv, prper=phprper, prvel=phprvel)
            amp1=10.**( (ampgr1-100.+amax) /20.)
            self.ftanparam.nfout1_1=nfout1
            arr1_1=np.concatenate((per1, tvis1, grvel1, phV1, ampgr1, ftrig1, snr1, wdth1, amp1))
            self.ftanparam.arr1_1=arr1_1.reshape(9, per1.size)
        else:
            amp1=10.**( (ampgr1-100.+amax) /20.)
            self.ftanparam.nfout1_1=nfout1
            arr1_1=np.concatenate((per1, tvis1, grvel1, phgr1, ampgr1, ftrig1, snr1, wdth1, amp1))
            self.ftanparam.arr1_1=arr1_1.reshape(9, per1.size)
            self.ftanparam.ampo_1=amp
        self.ftanparam.ncol_1, self.ftanparam.nrow_1 = amp.shape
        self.ftanparam.tamp_1=tamp = (nb-1)*dt+tb
        return
    
    def _compare_arr(self, data1, data2):
        plt.plot(data1-data2, '-y')
        plt.plot(data1, '^r')
        plt.plot(data2, '*b')
        plt.show()
        return
    
    def _compare_arr2(self, d1, d2, number):
        data1=d1[:, number]
        data2=d2[:data1.size, number]
        plt.plot(data1-data2, '-y')
        plt.plot(data1, '^r')
        plt.plot(data2, '*b')
        plt.show()
        return
    
    def taper(self, nb, ne, ntapb, ntape):
        omb = np.pi/ntapb
        ome = np.pi/ntape
        ncorr = int(ne+ntape)
        npts=self.stats.npts
        if ncorr>npts:
            ncorr=npts
        dataTapered=np.append(self.data[:ncorr], np.zeros( npts-ncorr ) )
        ##################################
        #zerp padding and cosine tapering
        ##################################
        # left end of the signal
        if nb-ntapb-1 > 0:
            dataTapered[:nb-ntapb-1]=0.
        if nb>ntapb:
            k=np.arange(ntapb+1)+nb-ntapb
            rwinb=(np.cos(omb*(nb-k))+1.)/2.
            dataTapered[nb-ntapb-1:nb]=rwinb*dataTapered[nb-ntapb-1:nb]
            sums = 2.*np.sum(rwinb)
        else:
            k=np.arange(nb)
            rwinb=(np.cos(omb*(nb-k))+1.)/2.
            dataTapered[:nb]=rwinb*dataTapered[:nb]
            sums = 2.*np.sum(rwinb)
        # right end of the signal
        if ne+ntape<npts:
            k=np.arange(ntape+1)+ne
            rwine=(np.cos(ome*(ne-k))+1.)/2.
            dataTapered[ne-1:ne+ntape] = dataTapered[ne-1:ne+ntape]*rwine
        elif ne < npts:
            k=np.arange(npts-ne+1)+ne
            rwine=(np.cos(ome*(ne-k))+1.)/2.
            dataTapered[ne-1:] = dataTapered[ne-1:]*rwine
        sums = sums+ne-nb-1
        c=np.sum(dataTapered[:ncorr])
        c=-c/sums
        # detrend
        if nb>ntapb:
            dataTapered[nb-ntapb-1:nb]=rwinb*c+dataTapered[nb-ntapb-1:nb]
        if ne+ntape<npts:
            dataTapered[ne-1:ne+ntape] = dataTapered[ne-1:ne+ntape] + rwine*c
        elif ne < npts:
            dataTapered[ne-1:] = dataTapered[ne-1:] + rwine*c
        dataTapered[nb:ne-1]=dataTapered[nb:ne-1]+c
        return dataTapered, ncorr
    
    def _fmax(self, amp, pha, ind, om, piover4 ):
        """parabolic interpolation of signal amplitude and phase, finding phase derivative
        """
        dt=self.stats.delta
        ind_l=ind-1
        ind_r=ind+1
        dd=amp[ind_l]+amp[ind_r]-2.*amp[ind]
        dd[dd==0]=-9999
        t=(amp[ind_l]-amp[ind_r])/dd/2.0 # shift from real maxima to ind_localmax
        t[dd==-9999]=0.
        a1=pha[ind_l]
        a2=pha[ind]
        a3=pha[ind_r]
        k1 = (a2-a1-om*dt)/2./np.pi
        k1=np.round(k1)
        a2 = a2-2.*k1*np.pi
        k2 = (a3-a2-om*dt)/2./np.pi
        k2=np.round(k2)
        a3 = a3-2.*k2*np.pi
        dph=t*(a1+a3-2.*a2)+(a3-a1)/2. # phase derivative
        tm=t*t*(amp[ind_l]+amp[ind_r]-2.*amp[ind])/2.+t*(amp[ind_l]-amp[ind_r])/2.+amp[ind] # calibrate amplitude maxima
        ph=t*t*(a1+a3-2.*a2)/2.+t*(a3-a1)/2.+a2+np.pi*piover4/4. # calibrate phase
        return dph, tm, ph, t
    
    def _phtovel(self, per, U, pha, npr, prper, prvel):
        """Convert observed phase to phase velocity
        """
        dist=self.stats.sac.dist
        omegaArr=2.*np.pi/per
        T=dist/U
        sU=1./U
        spl=scipy.interpolate.CubicSpline(prper, prvel)
        Vpred=spl(per[-1])
        phpred = omegaArr[-1]*(T[-1]-dist/Vpred)
        k=round((phpred -pha[-1])/2.0/np.pi)
        phV=np.zeros(U.size)
        phV[-1] = dist/(T[-1]-(pha[-1]+2.*k*np.pi)/omegaArr[-1])
        n=omegaArr.size
        for i in xrange(n-1):
            m=n-i-2
            Vpred =1/(((sU[m]+sU[m+1])*(omegaArr[m]-omegaArr[m+1])/2.+omegaArr[m+1]/phV[m+1])/omegaArr[m])
            phpred = omegaArr[m]*(T[m]-dist/Vpred)
            k = round((phpred -pha[m])/2.0/np.pi)
            phV[m] = dist/(T[m]-(pha[m]+2.0*k*np.pi)/omegaArr[m])
        return phV
    
    def _pred_cur(self, om0):
        """create phase prediction curve by group velocity, will be used for phase matched filter
        """
        pred=self.ftanparam.arr2_1[1:3,:]
        dist=self.stats.sac.dist
        x=2*np.pi/pred[0,::-1]
        y=dist/pred[1,::-1]
        ind_x=np.argsort(x)
        x=x[ind_x]
        y=y[ind_x]
        spl=scipy.interpolate.CubicSpline(x, y)
        gt0=spl(om0)
        y=y-gt0
        spl=scipy.interpolate.CubicSpline(x, y)
        return gt0, spl
        
    def _tapers(self, omb, ome, dom, alpha, ns):
        """spectra tapering
        """        
        om2d=omb/dom
        tresh=0.5
        wd = max(16., om2d*np.sqrt(tresh/alpha) )
        om1 = int(round(max(1, om2d-wd/2)))
        om2 = int(round(min(ns*1, om1+wd)))
        ampdom=np.zeros(ns)
        iArr1=np.arange(float(om2-om1+1))+om1
        ampdom[om1-1:om2]=(1.-np.cos(np.pi/(om2-om1)*(iArr1-om1)))/2.
        om3d = ome/dom
        wd = max(16., om3d*np.sqrt(tresh/alpha))
        om4 = int(round(min(ns*1, om3d+wd/2)))
        om3  = int(round(max(1, om4-wd)))
        iArr2=np.arange(float(om4-om3+1))+om3
        iArr2=iArr2[::-1]
        ampdom[om3-1:om4]=(1.-np.cos(np.pi/(om4-om3)*(iArr2-om3)))/2.
        ampdom[om2-1:om3]=1.
        omdom=np.arange(ns)*dom
        omstart = omb
        inds = om1
        inde = om4
        return omstart, inds, inde, omdom, ampdom
    
    def _tgauss(self, fsnr, gt0, dw, n, fmatch, seis):
        """taper phase matched signal
        """        
        ss=seis.copy()
        dt=self.stats.delta
        t0=self.stats.sac.b
        nc=round(gt0/dt)+1
        smax=np.abs(seis)
        ism=smax.argmax()
        sm=smax[ism]
        local_le=argrelextrema(smax, np.less_equal)[0]
        local_e=argrelextrema(smax, np.equal)[0]
        ind_localminima=np.setxor1d(local_le, local_e, assume_unique=True)
        ind_left=ind_localminima[ind_localminima<ism]
        ind_right=ind_localminima[ind_localminima>ism]
        val_left=smax[ind_left]
        val_right=smax[ind_right]
        nnnl=0
        if ind_left.size!=0:
            temp_nnnl=ind_left[((ism-ind_left)*dt>5.)*(val_left<fsnr*sm)]
            if temp_nnnl.size!=0:
                nnnl=temp_nnnl[-1]
                if temp_nnnl.size > 1:
                    nnl=temp_nnnl[-2]
                else:
                    nnl=0
        nnnr=0
        if ind_right.size!=0:
            temp_nnnr=ind_right[((ind_right-ism)*dt>5.)*(val_right<fsnr*sm)]
            if temp_nnnr.size!=0:
                nnnr=temp_nnnr[0]
                if temp_nnnr.size > 1:
                    nnr=temp_nnnr[1]
                else:
                    nnr=n-1        
        if nnnr!=0 and nnnl!=0:
            nn = max(abs(ism-nnnl), abs(ism-nnnr))
            nnn = max(abs(nnnl-nnl), abs(nnnr-nnr))
            nnnl = ism -nn
            nnl = nnnl-nnn
            nnnr = ism +nn
            nnr = nnnr+nnn
        tresh = np.log(sm)-24.
        if nnnl!=0:
            nnl = int(round((nnl-ism)*fmatch))+ism
            nnnl = int(round((nnnl-ism)*fmatch))+ism
            nnl = max(0, nnl)
            nnnl = max(0, nnnl)
            freq =(nnnl-nnl)+1
            iArr=np.arange(nnnl+1.)
            tre=-(iArr-nnnl)/freq*(iArr-nnnl)/freq/2.
            temp_ss=ss[:nnnl+1]
            temp_ss[tre>tresh]=temp_ss[tre>tresh]*(np.exp(tre))[tre>tresh]
            temp_ss[tre<=tresh]=0+0j
            ss[:nnnl+1]=temp_ss
        if nnnr!=0:
            nnr  = int(round((nnr-ism)*fmatch))+ism
            nnnr = int(round((nnnr-ism)*fmatch))+ism
            nnr  = min(n-1, nnr)
            nnnr = min(n-1, nnnr)
            freq = (nnr-nnnr)+1
            iArr = np.arange(float(n-nnnr))+nnnr+1
            tre  = -(iArr-nnnr-1)/freq*(iArr-nnnr-1)/freq/2.
            temp_ss=ss[nnnr:]
            temp_ss[tre>tresh]=temp_ss[tre>tresh]*(np.exp(tre))[tre>tresh]
            temp_ss[tre<=tresh]=0+0j
            ss[nnnr:]=temp_ss

        return ss
        
    def aftanf77(self, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, tmax=30.0, tresh=20.0,
            ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0, nfin=64, npoints=3, perc=50., phvelname='', predV=np.array([])):
        """ (Automatic Frequency-Time ANalysis) aftan analysis:
        ===========================================================================================================
        Input Parameters:
        pmf        - flag for Phase-Matched-Filtered output (default: True)
        piover4    - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin       - minimal group velocity, km/s
        vmax       - maximal group velocity, km/s
        tmin       - minimal period, s
        tmax       - maximal period, s
        tresh      - treshold for jump detection, usualy = 10, need modifications
        ffact      - factor to automatic filter parameter, usualy =1
        taperl     - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        snr        - phase match :q
        parameter, spectra ratio to determine cutting point for phase matched filter
        fmatch     - factor to length of phase matching window
        nfin       - number of initial period points
        npoints    - number of continuous points in jump correction
        perc       - output segment
        phvelname  - predicted phase velocity file name
        predV      - predicted phase velocity curve, period = predV[:, 0],  Vph = predV[:, 1]
        
        Output:
        self.ftanparam, a object of ftanParam class, to store output aftan results
        ===========================================================================================================
        """
        if not isaftanf77:
            raise AttributeError('fortran77 aftan not imported correctly!')
        # preparing for data
        try:
            self.ftanparam
        except:
            self.init_ftanParam()
        try:
            dist=self.stats.sac.dist
        except:
            dist, az, baz=obspy.geodetics.base.gps2dist_azimuth(self.stats.sac.evla, self.stats.sac.evlo,
                                self.stats.sac.stla, self.stats.sac.stlo) # distance is in m
            self.stats.sac.dist=dist/1000.
            dist=dist/1000.
        nprpv = 0
        phprper=np.zeros(300)
        phprvel=np.zeros(300)
        if predV.size != 0:
            phprper=predV[:,0]
            phprvel=predV[:,1]
            nprpv = predV[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        elif os.path.isfile(phvelname):
            # print 'Using prefile:',phvelname
            php=np.loadtxt(phvelname)
            phprper=php[:,0]
            phprvel=php[:,1]
            nprpv = php[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        else:
            warnings.warn('No predicted dispersion curve for:'+self.stats.network+'.'+self.stats.station, UserWarning, stacklevel=1)
        tempsac=self.copy()
        tb=self.stats.sac.b
        length=len(tempsac.data)
        if length>32768:
            warnings.warn('Length of seismogram is larger than 32768!', UserWarning, stacklevel=1)
            nsam=32768
            tempsac.data=tempsac.data[:nsam]
            tempsac.stats.e=(nsam-1)*tempsac.stats.delta+tb
            sig=tempsac.data
        else:
            sig=np.append(tempsac.data, np.zeros( 32768-tempsac.data.size, dtype='float64' ) )
            nsam=int( float (tempsac.stats.npts) )### for unknown reasons, this has to be done, nsam=int(tempsac.stats.npts)  won't work as an input for aftan
        dt=tempsac.stats.delta
        # Start to do aftan utilizing fortran 77 aftan
        self.ftanparam.nfout1_1,self.ftanparam.arr1_1,self.ftanparam.nfout2_1,self.ftanparam.arr2_1,self.ftanparam.tamp_1, \
                self.ftanparam.nrow_1,self.ftanparam.ncol_1,self.ftanparam.ampo_1, self.ftanparam.ierr_1= aftan.aftanpg(piover4, nsam, \
                    sig, tb, dt, dist, vmin, vmax, tmin, tmax, tresh, ffact, perc, npoints, taperl, nfin, snr, nprpv, phprper, phprvel)
        return

    def plot_group(self, title=''):
        """
        Plot ftan diagram: plot the group velocity dispersion curve
        ====================================================================
        title - title of the figure
        ====================================================================
        """
        try:
            fparam=self.ftanparam
            if fparam.nfout1_1==0:
                return "Error: No Basic FTAN parameters!"
            dt=self.stats.delta
            dist=self.stats.sac.dist
            v1=dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
            ampo_1=fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
            obper1_1=fparam.arr1_1[1,:fparam.nfout1_1]
            gvel1_1=fparam.arr1_1[2,:fparam.nfout1_1]
            phvel1_1=fparam.arr1_1[3,:fparam.nfout1_1]
            plb.figure()
            ax = plt.subplot()
            # p=plt.pcolormesh(obper1_1, v1, ampo_1, cmap='gist_rainbow',shading='gouraud')
            ax.plot(obper1_1, gvel1_1, '*r') #
            #cb = plt.colorbar(p, ax=ax)
            Tmin1=obper1_1[0]
            Tmax1=obper1_1[fparam.nfout1_1-1]
            vmin1= v1[fparam.ncol_1-1]
            vmax1=v1[0]
            plt.axis([Tmin1, Tmax1, vmin1, vmax1])
            plt.xlabel('Period(s)')
            plt.ylabel('Velocity(km/s)')
            plt.title(title)
        except AttributeError:
            print 'Error: FTAN Parameters are not available!'
        return
    
    def get_snr(self, ffact=1.):
        fparam=self.ftanparam
        dist=self.stats.sac.dist
        begT=self.stats.sac.b
        endT=self.stats.sac.e
        dt=self.stats.delta
        if fparam.nfout2_2!=0:
            o_per=fparam.arr2_2[1,:]
            g_vel=fparam.arr2_2[2,:]
            snrArr=np.ones(o_per.size)*-1.
            for i in xrange(fparam.nfout2_2):
                if g_vel[i]<0 or o_per[i]<0: continue
                filtered_tr=self.gaussian_filter_aftan(1./o_per[i], ffact=ffact)
                minT = dist/g_vel[i]-o_per[i]/2.
                maxT = dist/g_vel[i]+o_per[i]/2.
                if(minT<begT): minT=begT
                if(maxT>endT): maxT=endT
                # Noise window
                minT = maxT + o_per[i] * 5. + 500.
                skipflag=False
                if( (endT - minT) < 50. ): skipflag=True
                elif( (endT - minT) < 1100. ):
                    maxT = endT - 10.
                else:
                    minT = endT - 1100.
                    maxT = endT - 100.
                if not skipflag:
                    ib = (int)((minT-begT)/dt)
                    ie = (int)((maxT-begT)/dt)+2
                    tempnoise=filtered_tr[ib:ie]
                    noiserms=np.sqrt(( np.sum(tempnoise**2))/(ie-ib-1.) )
                    amp=self.ftanparam.arr2_2[7,i]
                    if noiserms!=0.: snrArr[i]=amp/noiserms
                    else: snrArr[i]=1.
            self.ftanparam.arr2_2=np.append(fparam.arr2_2, snrArr)
            self.ftanparam.arr2_2=self.ftanparam.arr2_2.reshape(9, o_per.size)
        return 
                
                    
        
    def gaussian_filter_snr(self, fcenter, fhlen=0.008):
        """
        Gaussian filter designed for SNR analysis, utilize pyfftw to do fft
        exp( (-0.5/fhlen^2)*(f-fcenter)^2 )
        ====================================================================
        Input parameters:
        fcenter - central period
        fhlen   - half length of Gaussian width
        ====================================================================
        """
        npts=self.stats.npts
        ns=1<<(npts-1).bit_length() # get an integer which is 2^(int)
        df=1.0/self.stats.delta/ns
        nhalf=ns/2+1
        fmax=(nhalf-1)*df
        if fcenter>fmax:
            fcenter=fmax
        alpha = -0.5/(fhlen*fhlen)
        F=np.arange(ns)*df
        gauamp = F - fcenter
        sf=np.exp(alpha*gauamp**2)
        if useFFTW:
            sp=pyfftw.interfaces.numpy_fft.fft(self.data, ns)
        else:
            sp=np.fft.fft(self.data, ns)
        filtered_sp=sf*sp
        filtered_sp[ns/2:]=0
        filtered_sp[0]/=2
        filtered_sp[ns/2-1]=filtered_sp[ns/2-1].real+0.j
        if useFFTW:
            filtered_seis=pyfftw.interfaces.numpy_fft.ifft(filtered_sp, ns)
        else:
            filtered_seis=np.fft.ifft(filtered_sp, ns)
        filtered_seis=2.*filtered_seis[:npts].real
        return filtered_seis
    
    def gaussian_filter_aftan(self, fcenter, ffact=1.):
        """
        Gaussian filter designed for SNR analysis, utilize pyfftw to do fft
        exp( (-0.5/fhlen^2)*(f-fcenter)^2 )
        ====================================================================
        Input parameters:
        fcenter - central period
        ffact   - factor to automatic filter parameter, usualy =1
        ====================================================================
        """
        npts=self.stats.npts
        ns=1<<(npts-1).bit_length()
        df=1.0/self.stats.delta/ns
        nhalf=ns/2+1
        fmax=(nhalf-1)*df
        alpha=ffact*20.
        if fcenter>fmax:
            fcenter=fmax
        omega0=2.*np.pi*fcenter
        omsArr=2.*np.pi*np.arange(ns)*df
        if useFFTW:
            sp=pyfftw.interfaces.numpy_fft.fft(self.data, ns)
        else:
            sp=np.fft.fft(self.data, ns)
        filtered_sp=_aftan_gaussian_filter(alpha=alpha, omega0=omega0, ns=ns, indata=sp, omsArr=omsArr)
        if useFFTW:
            filtered_seis=pyfftw.interfaces.numpy_fft.ifft(filtered_sp, ns)
        else:
            filtered_seis=np.fft.ifft(filtered_sp, ns)
        filtered_seis=2.*filtered_seis[:npts].real
        return filtered_seis
    
    
    
