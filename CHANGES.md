### v0.4.5 — 8 December 2021
- **Breaking change:** `bruges.attribute.similarity` has a new interface and computes its result(s) differently. I added several seismic coherence attributes as options to `bruges.attribute.similarity`. Methods include Marfurt coherence, that of Gersztenkorn, A., and K. J. Marfurt, 1999, and gradient structure tensor coherence (T Randen et al., 2000, 70th Annual International Meeting, SEG, Expanded Abstracts, 668-671.). Choose the method with the 'kind' argument. These functions were written by Joe Kington and published in The Leading Edge in the June 2015 issue.
- Complex trace attributes `instantaneous_amplitude` (envelope), `instantaneous_phase` and `instantaneous_frequency`, along with the `quadrature`, were added to `bruges.attribute`. 
- Fixed a bug in `bruges.models.wedge` that led to some wedge models not having the correct total size in the 'depth' dimension (e.g. sometimes being off by one pixel.)
- Now tested on Python 3.8, 3.9 (should have been working, just catching up!). Removing Python 3.4 and 3.5 from testing.
- Coming soon: dynamic-time-warp-assisted panel interpolation.

### v0.4.4 — 27 September 2021
- `bruges.reflection.acoustic_reflectivity()` and `bruges.reflection.reflectivity()` now allow you to pass the `axis` to compute on (axis 0 by default, which results in a reflectivity series for each _column_ in a 2D array). You can also set a 'mode', which is currently 'valid' and results in one less sample in the `axis` dimension. In the next release, v0.5, it will be 'same' by default, so the output is the same shape as the input.

### v0.4.3 — 30 June 2021
- New in the `models` submodule: `panel` allows you to pass in two or more 1D arrays, which can be of different lengths, and linearly interpolate in between them, sample by sample. Note: there is no attempt to correlate them; the assumption is that the top and bottom correlate and everything in between correlates in a linear way.
- New in `models.wedge`: varying net:gross along the 'breadth' (3rd) dimension of a wedge model. The wedge must be 'binary'; that is, the middle (wedge) layer must only contain 2 components. The 3rd dimension will then show pinch-out geometries that vary from one lithology to the other.
- It is now possible to choose `power` and `root` shapes for `models.wedge`, as well as `linear` and `sigmoid`.
- Fixed a bug in `models.wedge` that caused some wedges to fail with single values for the `width` argument.

### v0.4.2 — 1 March 2021
- The Ormsby wavelet now has the option of passing in the relative power (in dB) of the f2 and f3 corner frequencies, e.g. `P=(0, -5)`. Default: `(0, 0)` (the conventional trapezoidal Ormsby bandpass filter). 
- Wavelets now have the option of returning an odd number of samples for 'even' time periods, like 0.128 s at 0.004 s sample interval. This used to return 32 samples; now if you set `sym` to `True`, you'll get 33 samples. **Future change:** From v0.5, this will be the default.
- Wavelet time shoud no longer suffer from floating point imprecision (previously common near t = 0). 
- You can now optionally pass a time series to a wavelet function, and the wavelet will be evaluated at those times. This overrides `duration` and `dt`, so you should pass `None` to those parameters (or get a warning).
- Added a new module, `models`, which contains a new function, `wedge`. For now this generates 2D wedge models; in the future it will also provide 3D models.
- **Future change:** Reminder: wavelets in v0.5 will no longer have `return_t=True` by default.

### v0.4.1 — 21 February 2021
- Moved `phase_rotate()` to `bruges.filters` and it now should handle 2D and 3D seismic correctly. You can also pass an array-like of phases to get a 'phase bank' of rotations, similar to how this works with frequencies and wavelet banks. 

### v0.4 — 27 November 2020
- **Breaking change:** fixed numerous minor issues with `attribute.energy()`, see [issue 78](https://github.com/agile-geoscience/bruges/issues/78). Note that it now expects time to be in the last dimension of the array. The function now runs on n-D data and is also about 15 times faster.
- Multiple fixes to the documentation, thanks especially to Jesper Dramsch and Adriana Gordon.
- Added the `filters.berlage()` wavelet, a causal, minimum phase wavelet good for marine airgun sources.
- Added the `filters.generalized()` wavelet, of which the Ricker is a special case, as defined by Wang 2015, Geophys J Int 203, p1172ff.
- Added a FutureWarning that wavelets in the next release will have `return_t=True` by default. In a later release, the option will likely disappear completely.
- Added `fluids.rho_water()`, `fluids.rho_brine()`, `fluids.v_water()`, and `fluids.v_brine()`, implementing equations from Batzle & Wang (1992).
- Added `filters.convolve()`, which implements multi-dimensional convolution for 1D and 2D reflectivities and wavelets (or wavelet banks), and for 3D reflectivity with 1D wavelets. `convolve(rc_panel, wavelet_bank)` results in a 3D synthetic, with axes of frequency, offset, and time.
- Added `util.convolve_many()` to implement a fast way to convolve a 1D kernel over multiple signals whose domain (time or distance) is in the last axis.
- The reflectivity equations now accept a tuple for `theta1`, as well as a scalar value (interpreted as one angle, in degrees), or a list or array (interpreted as an array of angles). If you pass a 2-tuple, it is interpreted as a start and stop for a linear space; the step size will be 1.0 degrees. A 3-tuple is interpreted as a start, stop, and step size (note, not number of steps). The intervals specified will be open; the end-points will be included. So `theta1=(0, 5)` is equivalent to `theta1=[0,1,2,3,4,5]` and `theta1=(0, 10, 2)` is the same as `theta1=[0,2,4,6]`.

### v0.3.4 — 13 October 2018
- Added the NMO equation from Leo Uieda's TLE tutorial. This is a work in progress.
- Implemented `generic` (a windowed function), `cosine`, and `gabor` wavelets in `filters`. The `cosine` and `gabor` filters are implemented via `generic`.
- Fixed a bug in `filters.sinc` that caused the wavelet to have the wrong amplitude with no taper. Implemented `sinc` using `generic`.
- Fixed a bug in `reflection.critical_angles` that was preventing valid values from being computed.
- Fixed a bug in `CoordTransform` that broke Python 2.7 compatibility.
- Fixed up the documentation.  

### v0.3.3 — 7 February 2018 
- Fixed some bugs in v0.3.2, including a problem with non-integer indices in the windowed mean that `backus()` uses.
- Improved the time and depth conversion functions by adding the ability to return the new basis array (i.e. the new time or depth vector).

### v0.3.2 — 7 February 2018 
- Fixed a bug in `bruges.filters` that was returning the results from integer arrays as integers, giving incorrect results in those cases. Fixed tests.
- Reflectivity equations in `reflection` module now work on vectors, so you can use ndarrays for both the Vp, Vs, and rho values, and the theta values. They are about 10 times faster than running a loop over elements; the Zoeppritz is over 100 times faster. The results put the angles in the first dimension, so you can simply index in to get one offset.
- The Zoeppritz solutions and the Aki–Richards approximations now return the complex reflectivity and therefore show post-critical behaviour.
- New reflection coefficient series function, `reflection.reflectivity()` makes it easier to make offset reflectivities from logs.
- New acoustic reflection coefficient series function, `acoustic_reflectivity()`.
- Added `critical_angles()` and `reflection_phase()` functions to make it easier to compute the PP and PS critical angle(s), and to get the phase of a post-critical reflection.
- Improvements to `reflection` module docs.
- Deprecating  `moving_avg_conv` and `moving_avg_fft`, and for now the function `moving_average()` will return the convolutional solution.
- Several more simple linear and non-linear filters in `bruges.filters`, including `median` (good for seismic horizons) and `mode` (good for waveform classification).
- You can no longer import 'all' from a module. This is a bad idea anyway, so I'm not calling it a breaking change.
- The wavelets `ricker()` and `sweep()` now return transposed matrices if you ask for a wavelet bank by providing several frequencies. This is so the wavelets are in the first dimension, so you get get one by simply indexing.
- The `ormsby()` wavelet now also works for a sequence of frequency tuples, returning a wavelet bank.
- Fixed a bug in the `sweep()` wavelet that caused a time shift in the sweep function. Also added the taper option to the `sweep()` wavelet.
- Added a `sinc()` wavelet, with a taper option to attenuate the sidelobes.
- Added `inverse_gardner`, and other density and velocity transforms, to `petrophysics`.
- Added `transform.v_rms()` (RMS velocity), `transform.v_avg()` (average velocity) and `transform.v_bac()` (naïve Backus average). These all operate in a 'cumulative' average-down-to sense, per Liner (2004), Elements of 3D Seismology.

### v0.3.1 — 21 January 2018
- Repaired the `rockphysics.elastic_impedance()` function, which had various issues, including not applying normalization correctly and not working on vectors.

### v0.3.0 — 12 January 2018
- **Breaking change:** `anisotropy.backus` now always returns `Vp0`, `Vs0`, and `rho`.
- Non-compatible change: `filters.wavelets.ricker` no longer normalizes the output, to avoid a problem with amplitude dependence on filter length.
- Functions returning tuples now return a `collections.namedtuple` so it's more obvious what you're getting back.
- New nonlinear filters (plus tests): `filters.snn()`, `filters.kuwahara()` and `filters.conservative()`.
- `sweep()` optionally returns time basis, like the other wavelets.
- `reflection.shuey()` optionally returns intercept and gradient instead of the terms. Set `return_gradient=True`.
- No more `shuey2()`, `shuey3()`, `bortfeld2()`, `bortfeld3()`. Use `shuey()` and `bortfeld()` instead.
- Zoeppritz calculation is now fully vectorized and 75% faster. Also added the full `scattering_matrix()` function to the list of available functions (it was previously a private function).
- Better tests of the `reflection` algorithms (lower numerical tolerance).
- Tidied up code.

----

### v0.2.3 — 27 November 2017
- Coordinate transformation for (x,y) and (inline, line).
- Fixed documentation.

### v0.2.2 — 23 October 2017
- Fixes to elastic impedance.
- Returning rho from Backus average.
- Fixed Travis.

### v0.2.1 — 7 August 2015
- Fixed spectrogram bug.
- Fixed moduli warnings.

### v0.2.0 — 13 July 2015
- Anisotropy equations.

----

### v0.1.6 — 11 April 2015
- Float conversion in reflectivity.

### v0.1.5 — 11 April 2015
- Added keywords and changed defaults to spectral decomp.
- Added rock moduli.
- Fixed divide by zero in time to depth.

### v0.1.4 — 31 August 2014
- Added noise utility functions and rock physics calculators.

### v0.1.3 — 19 June 2014
- Made time to depth compatible for TWT.

### v0.1.2 — 13 June 2014
- Added data edge handling in time to depth conversion.

### v0.1.1 — December 12 2013
- Added AVO reflection models.
- Added time to depth conversions
- Added dip filters.
- Added spectral decomposition.
- Updated APIs to use time intervals and wavelet banks.

### v0.1.0 — 16 October 2013
- Initial release, a re-write of the agilegeo library.
