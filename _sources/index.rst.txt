:hide-toc:

.. container:: noclass
   :name: forkongithub

   `Fork on GitHub <https://github.com/agile-geoscience/bruges>`_

===========
bruges is a
===========

.. image:: http://agile.geosci.ai/bruges.png

In other words, it's just a load of functions that implement important equations in (mostly seismic) geophysics, from Aki-Richards to Zoeppritz.


Installation
------------

.. toctree::
    :caption: Quick start

To install ``bruges``, simply type the following into a terminal::

    pip install bruges

To create an Ormsby Wavelet::

    from bruges.filters.wavelets import ormbsy:
    w = ormsby(duration=0.5, dt=0.002, f=[5, 10, 40, 80])


User Guide
----------

    Coming soon!


API reference
-------------

.. toctree::
    :maxdepth: 3
    :caption: API reference

    bruges


Other resources
---------------

.. toctree::
    :maxdepth: 1
    :caption: Other resources

    development
    authors
    license
    changelog
    contributing


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
    :caption: Project links
    :hidden:

    PyPI releases <https://pypi.org/project/bruges/>
    Code in GitHub <https://github.com/agile-geoscience/bruges>
    Issue tracker <https://github.com/agile-geoscience/bruges/issues>
    Community guidelines <https://code.agilescientific.com/community>
    Agile's software <https://code.agilescientific.com>
    Agile's website <https://www.agilescientific.com>


.. include:: source/intro.rst

.. include:: source/status.rst

.. include:: source/install.rst
