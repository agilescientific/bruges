# bruges is a

![Bruges](http://agile.geosci.ai/bruges.png)

In other words, it's just a load of functions that implement important equations in (mostly seismic) geophysics, from Aki-Richards to Zoeppritz.

[![Run tests](https://github.com/agilescientific/bruges/actions/workflows/run-tests.yml/badge.svg)](https://github.com/agilescientific/bruges/actions/workflows/run-tests.yml)
[![Build docs](https://github.com/agilescientific/bruges/actions/workflows/build-docs.yml/badge.svg)](https://github.com/agilescientific/bruges/actions/workflows/build-docs.yml)
[![PyPI version](https://img.shields.io/pypi/v/bruges.svg)](https://pypi.python.org/pypi/bruges/)
[![PyPI versions](https://img.shields.io/pypi/pyversions/bruges.svg)](https://pypi.org/project/bruges//)
[![PyPI license](https://img.shields.io/pypi/l/bruges.svg)](https://pypi.org/project/bruges/)


## Quick start

Install with:

```shell
pip install bruges
```

Make a trapezoidal wavelet like:

```python
import bruges as bg
w, t = bg.filters.ormsby(duration=0.256, dt=0.002, f=[5, 10, 40, 80])
```

This produces two arrays: amplitude `w` and time `t`.


## Links

- [Documentation](https://code.agilescientific.com/bruges)
- [Issue Tracker](https://github.com/agilescientific/bruges/issues/)
- [PyPi](http://pypi.python.org/pypi/bruges/)
- [Agile's website](http://www.agilescientific.com)

![Bruges rooves](https://www.dropbox.com/s/tzvi22ujq6rozdb/bruges_long_rooves.png?raw=1)
