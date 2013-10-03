from traits.api import HasTraits, Int, Event
from traitsui.api import Item
from numpy import zeros

def integrate_spectra( data, b1, b2 ):
    """
    Integrates across each slice of a spectrogram
    """
    output = data[ :, b1:b2].sum( axis = 1 )
    
    return output
 
class SpectralDecomposition( HasTraits ): 
    
    b1_low = Int
    b1_high = Int
    
    b2_low = Int
    b2_high = Int

    b3_low = Int
    b3_high = Int
    
    updated = Event
    
    traits_view=( Item( 'b1_low'),
                  Item( 'b1_high'),
                  Item( 'b2_low'),
                  Item( 'b2_high'),
                  Item( 'b3_low'),
                  Item( 'b3_high') )
                  
    def __init__( self, b1_low=5, b1_high=10, b2_low=10, b2_high=20,
                        b3_low=20, b3_high=40):
        
        super( SpectralDecomposition, self ).__init__()
        
        self.b1_low = b1_low
        self.b1_high = b1_high
        self.b2_low = b2_low
        self.b2_high = b2_high
        self.b3_low = b3_low
        self.b3_high = b3_high
        
    def _on_anytrait_changed( self ):
        self.updated = True
      
      
    def compute( self, data ):
        
        output = zeros( ( data.shape[0], data.shape[1], 3 ) )
        
        for spectra in range( data.shape[-2] ):
            
            b1 = integrate_spectra( data[ :, spectra, :], self.b1_low, 
                                    self.b1_high )
            b2 = integrate_spectra( data[ :, spectra, :], self.b2_low, 
                                    self.b2_high )
            b3 = integrate_spectra( data[ :, spectra, :], self.b3_low, 
                                    self.b3_high )
                                    
            output[ :, spectra, 0 ] = b1
            output[ :, spectra, 1 ] = b2
            output[ :, spectra, 2 ] = b3
        return( output )
        
        