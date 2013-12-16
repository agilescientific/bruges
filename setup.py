from setuptools import setup

    
setup(
    name='agilegeo',
    version=open('version.txt').read().rstrip(),
    author='Agile Geoscience',
    author_email='hello@agilegeoscience.com',
    packages=['agilegeo',
              'agilegeo.wavelet',
              'agilegeo.attribute',
              'agilegeo.avo',
              'agilegeo.util'],
    description='Useful geophysics functions',
    long_description=open('README.rst').read(),
    url='http://pypi.python.org/pypi/agilegeo/',
    license='LICENSE.txt'
    )
