### v0.3.2 — February 2018 
- Fixed a bug in `bruges.filters` that was returning the results from integer arrays as integers, giving incorrect results in those cases. Fixed tests.
- Reflectivity equations in `reflection` module now work on vectors. They are about 10 times faster than running a loop over elements. The results put the angles in the first dimension, so you can simply index in to get one offset.
- The Zoeppritz solutions and the Aki–Richards approximations now return the complex reflectivity and therefore show post-critical behaviour.
- New reflection coefficient series function, `reflection.reflectivity()` makes it easier to make offset reflectivities from logs.
- New acoustic reflection coefficient series function, `acoustic_reflectivity()`.
- Improvements to `reflection` module docs. 
- Deprecating  `moving_avg_conv` and `moving_avg_fft`, and for now the function `moving_average()` will return the convolutional solution.
- Several more simple linear and non-linear filters in `bruges.filters`. They are n-dimensional where possible. (One day maybe I'll have a crack at n-D SNN and Kuwahara filters.)
- You can no longer import 'all' from a module. This is a bad idea anyway,
so I'm not calling it a breaking change.
- The wavelets `ricker()` and `sweep()` now return transposed matrices if you ask for a wavelet bank by providing several frequencies. This is so the wavelets are in the first dimension, so you get get one by simply indexing.
- Added `inverse_gardner`, and other density and velocity transforms, to `petrophysics`.
- Added `transform.v_rms()` (RMS velocity), `transform.v_avg()` (average velocity) and `transform.v_bac()` (naïve Backus average). These all operate in a 'cumulative' average-down-to sense, per Liner (2004), Elements of 3D Seismology.

### v0.3.1 — 21 January 2018
- Repaired the `rockphysics.elastic_impedance()` function, which had various issues, including not applying normalization correctly and not working on vectors.

### v0.3.0 — 12 January 2018
- Breaking change: `anisotropy.backus` now always returns `Vp0`, `Vs0`, and `rho`.
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
