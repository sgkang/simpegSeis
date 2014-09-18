.. image:: simpeg-logo.png
   :width: 300 px
   :alt: SimPEG
   :align: center

SimPEG (Simulation and Parameter Estimation in Geophysics) is a python package for simulation and gradient based parameter estimation in the
context of geoscience applications. simpegSeis uses SimPEG as the framework for the forward modeling and inversion of seismic geophysical problem.

We consider two fundamental governing equations for seismic wave propagation: acoustic and elastic wave equations in both time and frequency domains. We initially discretize those problems to compute forward responses. Obvious next step is to solve inverse problems

A fundamental purpose of geophysical application is to better image the subsurface of the earth. One geophysical method only provide us a specific information about related material properties to that method. Therefore, to successfully achieve this goal, we should integrate different geophyisical applications. We believe that we, who have different geophysical backgrounds (Sesmic, EM, Magnetics, Gravity, ...), should cooperate. Therefore, we welcome that any people who want to contribute this project, simpegSeis. 

.. raw:: html
    :file: examples/refraction.html

Acoustic wave
=============

.. toctree::
   :maxdepth: 2

   api_Acous

Elastic wave
============


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

