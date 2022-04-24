"""
GLCM (gray level co-occurrence matrices)
:copyright: 2022 Software Underground 
:license: Apache 2.0
"""

from bruges.util.util import patchify, unpatchify
from skimage.feature.texture import greycomatrix, greycoprops
from skimage.transform import resize, warp
import skimage
from scipy.ndimage import convolve
import numpy as np
from tqdm import tqdm

def main():
    pass

def glcm2d(img, vmin=None, vmax=None, levels=8, kernel_size=5, distance=1.0, angle=0.0, feature='contrast', method='skimage'):
    '''
    Compute statistical properties from GLCM (gray level co-occurrence matrices) on seismic data.
    https://subsurfwiki.org/wiki/Haralick_texture
    Algorithm :
        1. Read data from an image, seismic horizon, or grid
        2. Compute GLCMs in a moving window
        3. Compute statistics of GLCMs, and put in new grid
        4. Result: various attributes, such as entropy, energy, contrast
    
    Source:
    GLCM
        Haralick, R (1979). Statistical and structural approaches to texture. Proceedings of the IEEE 67, 786–804.
    GLCM using skimage
        Stéfan van der Walt, Johannes L. Schönberger, Juan Nunez-Iglesias, François Boulogne, Joshua D. Warner, Neil Yager, Emmanuelle Gouillart, Tony Yu and the scikit-image contributors. scikit-image: Image processing in Python. PeerJ 2:e453 (2014) https://doi.org/10.7717/peerj.453
        https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_glcm.html?highlight=glcm
    GLCM on seismic data
        Marcilio Castro de Matos, Malleswar (Moe) Yenugu, Sipuikinene Miguel Angelo, and Kurt J. Marfurt, (2011), "Integrated seismic texture segmentation and cluster analysis applied to channel delineation and chert reservoir characterization," GEOPHYSICS 76: P11-P21.https://doi.org/10.1190/geo2010-0150.1
        https://library.seg.org/doi/abs/10.1190/geo2010-0150.1?journalCode=gpysa7#
    Fast GLCM 
        tzm030329, Fast GLCM, (2018), GitHub repository, https://github.com/tzm030329/GLCM
    
    Args:
    ----------
    img: 2darray, shape=(h,w), dtype=np.uint8 / float
        input image
    vmin: float
        minimum value of input image
    vmax: float
        maximum value of input image
    levels: int
        number of grey-levels of GLCM
    kernel_size: int
        Patch size to calculate GLCM around the target pixel
    distance: float
        pixel distance offsets (1.0, 2.0, and etc.)
    angle: float
        pixel angles (0.0, 30.0, 45.0, 90.0, and etc.)
    feature: string
        statistical feature from GLCM ('contrast','dissimilarity','homogeneity','ASM','energy')
    method: string
        method to use : ('skimage', 'fast')
            'skimage': default,
            'fast': Fastest GLCM calculation using opencv warp affine and 2D Convolution
        
    Returns:
    -------
    result: 2darray
        2D Grey-level co-occurrence matrix statistical feature array the same shape as the input array.

    '''
    if (vmin == None) or (vmax == None):
        mi, ma = img.min(), img.max()
    else:
        mi, ma = vmin, vmax
    ks = kernel_size
    h,w = img.shape

    # digitize
    bins = np.linspace(mi, ma+1, levels+1)
    gl1 = np.digitize(img, bins) - 1


    if method == 'skimage':        
        patches = patchify(gl1, (kernel_size, kernel_size), step=1)
        glcmstat = np.zeros_like(patches[:,:,0,0])
        with tqdm(total=patches.shape[0]*patches.shape[1], desc="GLCM statistic calculation") as pbar:
            for i in range(patches.shape[0]):
                for j in range(patches.shape[1]):
                    glcm2 = greycomatrix(patches[i,j,:,:], distances=[int(distance)], angles=[int(angle)], levels=levels,
                                            #symmetric=True, 
                                            normed=True
                                            )
                    result2 = greycoprops(glcm2, feature)
                    glcmstat[i,j] = result2[0][0]
                    pbar.update(1)
        pbar.close()        
        result = resize(glcmstat,(gl1.shape[0], gl1.shape[1]))
        
    if method == 'fast':
        # make shifted image
        dx = distance*np.cos(np.deg2rad(angle))
        dy = distance*np.sin(np.deg2rad(-angle))
        #mat = np.array([[1.0,0.0,-dx], [0.0,1.0,-dy]], dtype=np.float32)
        mat2 = np.array([[1.0,0.0,-dx], [0.0,1.0,-dy], [0.0,0.0,1.0]], dtype=np.float32)
        gl2 = warp(gl1.astype(float), mat2, mode='edge')
        # make glcm
        glcm = np.zeros((levels, levels, h, w), dtype=np.uint8)
        for i in range(levels):
            for j in range(levels):
                mask = ((gl1==i) & (gl2==j))
                glcm[i,j, mask] = 1

        kernel = np.ones((ks, ks), dtype=np.uint8)
        for i in range(levels):
            for j in range(levels):
                glcm[i,j] = convolve(glcm[i,j], kernel, mode='nearest')

        glcm = glcm.astype(np.float32)

        #select feature
        if feature == 'contrast':
            result = np.zeros((h,w), dtype=np.float32)
            for i in range(levels):
                for j in range(levels):
                    result += glcm[i,j] * (i-j)**2 
                    
        if feature == 'dissimilarity':
            result = np.zeros((h,w), dtype=np.float32)
            for i in range(levels):
                for j in range(levels):
                    result += glcm[i,j] * np.abs(i-j)
                    
        if feature == 'homogeneity':    
            result = np.zeros((h,w), dtype=np.float32)
            for i in range(levels):
                for j in range(levels):
                    result += glcm[i,j] / (1.+(i-j)**2)
                    
        if feature == 'ASM':
            result = np.zeros((h,w), dtype=np.float32)
            for i in range(levels):
                for j in range(levels):
                    result  += glcm[i,j]**2
                    
        if feature == 'energy':
            result = np.zeros((h,w), dtype=np.float32)
            for i in range(levels):
                for j in range(levels):
                    result  += glcm[i,j]**2
            result = np.sqrt(result)
    #normalizing        
    result =  result/(result.max()/1.0)        
    return result
    
    
if __name__ == '__main__':
    main()

    img = skimage.data.camera()
    result = glcm2d(data, levels=8, kernel_size=3, distance=1.0, angle=0.0, feature='dissimilarity', method='fast')
