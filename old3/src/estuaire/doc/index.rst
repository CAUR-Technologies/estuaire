.. Delta documentation master file, created by
   sphinx-quickstart on Mon Sep 27 16:09:46 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Estuaire's documentation!
====================================

Contents:

.. toctree::
   :maxdepth: 2

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

License
========
All the libraries/source code and documentation are property of Jean-Pascal Mercier and Advanced GeoScience Imaging Solution Ltd.


Preface
========
This document is in a preliminary state.

About
======



Package Contents
================

Installation
=============

This document contains the installation instructions for estuaires Seismic Flow Control System, Eikonal Solver and Raytracing software. The whole system contains 3 interdependent modules. Two of them entirely written in Python, the other mixes C++ and Python for performance boost in non-linear algorithms.


Requirements
-------------
Since all the modules main API are in Python, you must obviously have Python installed on your system.

The Raytracer and Eikonal Solver C++/Python bridge is written in Cython and make good use of Mako template engine. Finally, the bridge also use Numpy/Scipy numerical libraries and development files.

The flow control system is built over scons.

Recap.:
    1) Python Development Files
    2) Cython <http://www.cython.org/>
    3) Mako Template Engine <http://www.makotemplates.org>
    4) Numpy/Scipy <http://scipy.org>
    5) SCons <http://www.scons.org>
    6) Argparse
    7) PyTables



.. note:: Currently, the only tested platform is GNU/Linux Operating System. This should probably work under any Unix flavors with little or no work. There is no reason why it could not be compiled/installed on Windows, althought it has never been tried or tested in any ways.

Eikonal Solver
---------------
This modules contains the low level C++ Eikonal Solver as well as the Frechet derivatives and Raytracing capabilities. Both the low level C++ implementation and the python binding are provided.

Eikonal-ng
==========
.. include:: /raid/jee/agsis/development/eikonal-ng/README.rst

SCons Subsytem
==============

One of the important feature of every scientific experiment is the possibility to reproduce the results. As a scientific software, Estuaire want make strong commitment towards producing scientifically sound computational results that can be replicate again and again.

At low level, Estuaire rely strongly on a powerful SCons subsystem.

SCons can be think about as a boosted GNU Make dependency handling system coupled with the strong, fast and extensible Python scripting language. This enables a complete flexibility and power to accomplish more.

Even though SCons was primarly developped as a building system, the basic dependency handling can be extended in a manner to accomplish any tasks stream that can be represented as a Directed Acyclic Graph (DAG).

Usage
=====

The site_scons environmenent is provided with basically 2 levels of API.


.. warning::
    .. include:: ../README_SCONS_2.rst

Lotic
------
Add a Simple Example of a SConstruct

.. automodule:: site_tools.Lotic


Slopes
-------
.. automodule:: slopes

Database Filter
----------------
.. automodule:: dbfilter
    :members: __doc__, CuboidFilter, DateFilter, StationFilter, EventFilter

