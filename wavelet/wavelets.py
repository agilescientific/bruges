# -*- coding: utf-8 -*-
import numpy as np
from enthought.traits.api import HasTraits,Array,Float, Range, Int, Event
from scipy.signal import hilbert
from enthought.traits.ui.api import View,Item, RangeEditor

def ricker( duration, dt, f ):
    """
    Also known as the mexican hat wavelet, models the function:
    A =  (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}

    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001, 0.002, 0.004)
    :params f: Center frequency of the wavelet (in Hz)

    :returns: The ricker wavelet with center frequency f sampled at t.
    """
    t = np.arange( -duration/2, duration/2 , dt) 
    pi2 = (np.pi ** 2.0)
    fsqr = f ** 2
    tsqr = t ** 2
    pft = pi2 * fsqr * tsqr
    A = (1 - (2 * pft)) * np.exp(-pft)

    return A
    
    
def sweep( duration, dt, f1, f2, method = 'linear', phi = 0, vertex_zero = True, 
    autocorrelate = True ):
    """
    Generates a linear frequency modulated wavelet (sweep)
    Does a wrapping of scipy.signal.chirp
    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001, 0.002, 0.004)
    :param f1: The start frequency of the wavelet.
    :param f2: The end frequency of the wavelet (in Hz)
    :keyword method: {'linear','quadratic','logarithmic'}, optional
    :keyword phi: float, phase offset in degrees
    :keyword vertex_zero: bool, optional
        This parameter is only used when method is ‘quadratic’. 
        It determines whether the vertex of the parabola that 
        is the graph of the frequency is at t=0 or t=t1.

    :returns: An LFM waveform.
    """
    
    t = np.arange( -duration/2, duration/2 , dt) 
    t0 = -duration/2
    t1 = duration/2
    from scipy.signal import chirp
    A = chirp( t , f1, t1, f2, method, phi, vertex_zero  )
    
    if autocorrelate:
        A = np.correlate(A, A, mode='same')
    
    return A


def ormsby(duration, dt, f1, f2, f3, f4):
    """
    The Ormsby wavelet requires four frequencies:
    f1 = low-cut frequency
    f2 = low-pass frequency
    f3 = high-pass frequency
    f4 = hi-cut frequency
    Together, the frequencies define a trapezoid shape in the spectrum.
    The Ormsby wavelet has several sidelobes, unlike Ricker wavelets
    which only have two, one either side.

    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001, 0.002, 0.004)
    :params f1: Low cut frequency
    :params f2: Low pass frequency
    :params f3: High pass frequency
    :params f4: High cut frequency

    :returns: A vector containing the ormsby wavelet
    """
    
    def numerator(f,t):
        return (np.sinc( f * t)**2) * ((np.pi * f) ** 2)

    pf43 = ( np.pi * f4 ) - ( np.pi * f3 )
    pf21 = ( np.pi * f2 ) - ( np.pi * f1 )
    
    t = np.arange( -duration/2, duration/2 , dt) 
    
    A = ( ( numerator(f4,t)/pf43 ) - ( numerator(f3,t)/pf43 ) -
         ( numerator(f2,t)/pf21 ) + ( numerator(f1,t)/pf21) )

    A /= np.amax( A )
    return A


def rotate_phase( w, phi ):
    """
    Performs a phase rotation of wavelet using:
    A = w(t)Cos(phi) + h(t)Sin(phi)
    Where w(t) is the wavelet and h(t) is it's hilbert transform.
    
    :params w: The wavelet vector.
    :params phi: The phase rotation angle (in Radians) to apply.

    :returns: The phase rotated signal.
    """

    # Get the analytic signal for the wavelet
    a = hilbert( w )

    A = (np.real( a ) * np.cos( phi ) +
         np.imag( a  ) * np.sin( phi ) )

    return A

class Wavelet( HasTraits ):
    """
    Abstract Class
    """
    wavelet = Array
    nsamps = Int
    duration = Float
    phase = Range( -180.,180.,0, mode = 'slider')
    updated = Event
    
    def __init__(self):
        super( Wavelet, self).__init__()
    
    

class Ricker( Wavelet ):
    
    nsamps = Int(512)
    duration = Float(0.2)
    center_frequency = Range( 1., 100. ,40)

   
    triats_view = View( Item( 'nsamps' ), 
                        Item( 'duration' ),
                        Item( 'center_frequency' ),
                        Item( 'phase' ),
                        buttons=['OK'] )
          
    def __init__(self):
        super( Ricker, self).__init__()
        self.updated = True
        self._compute()
                      
    def _compute( self ):
        
        wavelet = ricker_alg( self.duration, 
                              self.nsamps, 
                              self.center_frequency )
        
        phase_rad = self.phase * np.pi / 180.                  
        self.wavelet = rotate_phase( wavelet,
                                     phase_rad )
        
                                   
    def _anytrait_changed( self, name ):
        if (( name not in ['wavelet', 'updated'] ) ):
            self._compute()
            self.updated = True

class Ormsby( Wavelet ):
    
    max_f = Float( 100.0 )
    min_f = Float( 0.0 )
    f1 = Range(low=1., value=10.)
    f2 = Range(low=1., value=30.)
    f3 = Range(low=1.,value=50.)
    f4 = Range(low=1., value=75.)
    
    phase = Range(-180.,180.,0., mode='slider' )
    nsamps = Int( 512 )
    duration = Float(0.2)
  
    traits_view = View( Item('nsamps'), Item( 'duration' ),
                        Item( 'f1', editor=RangeEditor( mode='slider',
                                                        low_name = 'min_f',
                                                        high_name='f2' )), 
                        Item( 'f2', editor=RangeEditor( mode='slider', 
                                                        low_name='f1',
                                                        high_name='f3' )),
                        Item( 'f3',editor=RangeEditor( mode='slider', 
                                                        low_name='f2',
                                                        high_name='f4' )),
                        Item( 'f4', editor=RangeEditor( mode='slider', 
                                                        low_name='f3',
                                                        high_name='max_f')),
                        Item( 'phase' ),
                        buttons = ['OK'])
    def __init__( self ):
        super( Ormsby, self ).__init__()
        self.updated=True
        self._compute()
        
    
    def _compute( self ):
        wavelet = ormsby_alg( self.duration, self.nsamps,
                              self.f1, self.f2, self.f3,
                              self.f4 )
        phase_rad = self.phase * np.pi / 180.
        self.wavelet = rotate_phase( wavelet, phase_rad )
                        
    def _anytrait_changed( self, name ):
        if (( name not in ['wavelet', 'updated'] ) ):
            self._compute()
            self.updated = True
        
    
    



    
