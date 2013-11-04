from distutils.core import setup

setup(
    name='agilegeo',
    version='0.1.0',
    author='Ben Bougher',
    author_email='ben.bougher@gmail.com',
    packages=['agilegeo','agilegeo.wavelet', 'agilegeo.attribute',
              'agilegeo.attribute.test'],
    description='Useful geophysics functions',
    long_description=open('README.txt').read(),
    url='http://pypi.python.org/pypi/agilegeo/',
    license='LICENSE.txt'
    )
