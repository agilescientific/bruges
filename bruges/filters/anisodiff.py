"""
:copyright: Alistair Muldal
:license: Unknown, shared on StackOverflow and Pastebin

Reference:
P. Perona and J. Malik.
Scale-space and edge detection using ansotropic diffusion.
IEEE Transactions on Pattern Analysis and Machine Intelligence,
12(7):629-639, July 1990.
<http://www.cs.berkeley.edu/~malik/papers/MP-aniso.pdf>

Original MATLAB code by Peter Kovesi
School of Computer Science & Software Engineering
The University of Western Australia
pk @ csse uwa edu au
<http://www.csse.uwa.edu.au>

Translated to Python and optimised by Alistair Muldal
Department of Pharmacology
University of Oxford
<alistair.muldal@pharm.ox.ac.uk>

June 2000  original version.
March 2002 corrected diffusion eqn No 2.
July 2012 translated to Python
"""
import numpy as np
import warnings


def anisodiff(img, niter=1, kappa=50, gamma=0.1, step=(1., 1.), option=1):
    """
    Anisotropic diffusion.

    Usage:
    imgout = anisodiff(im, niter, kappa, gamma, option)

    Arguments:
            img    - input image
            niter  - number of iterations
            kappa  - conduction coefficient 20-100 ?
            gamma  - max value of .25 for stability
            step   - tuple, the distance between adjacent pixels in (y,x)
            option - 1 Perona Malik diffusion equation No 1
                     2 Perona Malik diffusion equation No 2

    Returns:
            imgout   - diffused image.

    kappa controls conduction as a function of gradient.  If kappa is low
    small intensity gradients are able to block conduction and hence
    diffusion across step edges.  A large value reduces the influence of
    intensity gradients on conduction.

    gamma controls speed of diffusion (you usually want it at a maximum of
    0.25)

    step is used to scale the gradients in case the spacing between
    adjacent pixels differs in the x and y axes

    Diffusion equation 1 favours high contrast edges over low contrast
    ones.
    Diffusion equation 2 favours wide regions over smaller ones.
    """

    # ...you could always diffuse each color channel independently if you
    # really want
    if img.ndim == 3:
        m = "Only grayscale images allowed, converting to 2D matrix"
        warnings.warn(m)
        img = img.mean(2)

    # initialize output array
    img = img.astype('float32')
    imgout = img.copy()

    # initialize some internal variables
    deltaS = np.zeros_like(imgout)
    deltaE = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    gS = np.ones_like(imgout)
    gE = gS.copy()

    for ii in xrange(niter):

        # calculate the diffs
        deltaS[:-1, :] = np.diff(imgout, axis=0)
        deltaE[:, :-1] = np.diff(imgout, axis=1)

        # conduction gradients (only need to compute one per dim!)
        if option == 1:
                gS = np.exp(-(deltaS/kappa)**2.)/step[0]
                gE = np.exp(-(deltaE/kappa)**2.)/step[1]
        elif option == 2:
                gS = 1./(1.+(deltaS/kappa)**2.)/step[0]
                gE = 1./(1.+(deltaE/kappa)**2.)/step[1]

        # update matrices
        E = gE*deltaE
        S = gS*deltaS

        # subtract a copy that has been shifted 'North/West' by one
        # pixel. don't as questions. just do it. trust me.
        NS[:] = S
        EW[:] = E
        NS[1:, :] -= S[:-1, :]
        EW[:, 1:] -= E[:, :-1]

        # update the image
        imgout += gamma*(NS+EW)

    return imgout


def anisodiff3(stack, niter=1, kappa=50, gamma=0.1, step=(1., 1., 1.), option=1):
    """
    3D Anisotropic diffusion.

    Usage:
    stackout = anisodiff(stack, niter, kappa, gamma, option)

    Arguments:
        stack  - input stack
        niter  - number of iterations
        kappa  - conduction coefficient 20-100 ?
        gamma  - max value of .25 for stability
        step   - tuple, the distance between adjacent pixels in (z,y,x)
        option - 1 Perona Malik diffusion equation No 1
                 2 Perona Malik diffusion equation No 2

    Returns:
        stackout   - diffused stack.

    kappa controls conduction as a function of gradient.  If kappa is low
    small intensity gradients are able to block conduction and hence
    diffusion across step edges. A large value reduces the influence of
    intensity gradients on conduction.

    gamma controls speed of diffusion (you usually want it at a maximum of
    0.25)

    step is used to scale the gradients in case the spacing between
    adjacent pixels differs in the x,y and/or z axes

    Diffusion equation 1 favours high contrast edges over low contrast
    ones.
    Diffusion equation 2 favours wide regions over smaller ones.
    """

    # ...you could always diffuse each color channel independently if you
    # really want
    if stack.ndim == 4:
        m = "Only grayscale stacks allowed, converting to 3D matrix"
        warnings.warn(m)
        stack = stack.mean(3)

    # initialize output array
    stack = stack.astype('float32')
    stackout = stack.copy()

    # initialize some internal variables
    deltaS = np.zeros_like(stackout)
    deltaE = deltaS.copy()
    deltaD = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    UD = deltaS.copy()
    gS = np.ones_like(stackout)
    gE = gS.copy()
    gD = gS.copy()

    for ii in range(niter):

        # calculate the diffs
        deltaD[:-1, :, :] = np.diff(stackout, axis=0)
        deltaS[:, :-1, :] = np.diff(stackout, axis=1)
        deltaE[:, :, :-1] = np.diff(stackout, axis=2)

        # conduction gradients (only need to compute one per dim!)
        if option == 1:
                gD = np.exp(-(deltaD/kappa)**2.)/step[0]
                gS = np.exp(-(deltaS/kappa)**2.)/step[1]
                gE = np.exp(-(deltaE/kappa)**2.)/step[2]
        elif option == 2:
                gD = 1./(1.+(deltaD/kappa)**2.)/step[0]
                gS = 1./(1.+(deltaS/kappa)**2.)/step[1]
                gE = 1./(1.+(deltaE/kappa)**2.)/step[2]

        # Update matrices.
        D = gD*deltaD
        E = gE*deltaE
        S = gS*deltaS

        # Subtract a copy that has been shifted 'Up/North/West' by one
        # pixel. Don't ask questions. Just do it. Trust me.
        UD[:] = D
        NS[:] = S
        EW[:] = E
        UD[1:, :, :] -= D[:-1, :, :]
        NS[:, 1:, :] -= S[:, :-1, :]
        EW[:, :, 1:] -= E[:, :, :-1]

        # update the image
        stackout += gamma*(UD+NS+EW)

    return stackout
