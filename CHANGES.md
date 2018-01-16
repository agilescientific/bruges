### v0.3.1 — January 2018
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
