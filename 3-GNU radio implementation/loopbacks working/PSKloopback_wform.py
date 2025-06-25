"""This script defines the pulse-shaping filter (root raised cosine) that is used in the psk modulation."""

import math
import numpy

#Nyquist pulse: without inter symbol interference
def nyq(Sps,ntaps):
    n=numpy.linspace(-int(ntaps/2), int(ntaps/2-1),ntaps) # symmetric time vector centered at 0
    h=numpy.sinc(n/Sps)  # scale sinc function
    return h

# Root Raised Cosine Filter: for transmit/receive matched filtering
# convolution with itself
def rrcos(Sps,ntaps,beta):
    if beta==0:
        h=nyq(Sps,ntaps) #ideal nyquist pulse
    else:
        h=ntaps*[0,]
        beta4=4.*beta
        for n in range(ntaps):
            k=n-ntaps/2. #  h[n] center tap at 0
            if k==0:
                h[n]=1+beta*(4./math.pi-1.)
            elif abs(k)==Sps/beta4:
                ha=(1.+2./math.pi)*math.sin(math.pi/beta4)
                hb=(1.-2./math.pi)*math.cos(math.pi/beta4)
                h[n]=(ha+hb)*beta/math.sqrt(2.)
            else:
                ks=k/Sps #rrc formula: time impulse response
                kspi=math.pi*ks
                Num=math.sin(kspi*(1-beta))+beta4*ks*math.cos(kspi*(1+beta))
                Den=kspi*(1.-(beta4*ks)**2)
                h[n]=Num/Den
    Amp=numpy.amax(h) #normalize by the max amplitude
    return h/Amp

#sinc filter (used for Nyquist filtering)
def sinc(Sps,ntaps):
    n=np.linspace(-int(ntaps/2), int(ntaps/2-1),ntaps)
    h=np.sinc(n/Sps) #sinc function
    return h

