.. _api_Acous:


.. math::

    \renewcommand{\div}{\nabla\cdot\,}
    \newcommand{\grad}{\vec \nabla}
    \newcommand{\curl}{{\vec \nabla}\times\,}
    \newcommand{\Acf}{{\mathbf A_c^f}}
    \newcommand{\Ace}{{\mathbf A_c^e}}
    \newcommand{\dcurl}{{\mathbf C}}
    \newcommand{\dgrad}{{\mathbf G}}
    \newcommand{\u}{\vec{u}}
    \renewcommand{\dudt}{\frac{\partial\u}{\partial t}}
    \renewcommand{\dphidt}{\frac{\partial\phi}{\partial t}}
    \renewcommand{\dsdt}{\frac{\partial s}{\partial t}}

Acoustic wave equation
**********************

Time domain acoustic wave equation in 1st order form:

.. math::

        \rho\frac{\partial \phi}{\partial t} = \nabla \cdot \vec{u}+\frac{\partial s}{\partial t}\delta(\vec{r}-\vec{r}_s)

        \mu^{-1}\frac{\partial \vec{u}}{\partial t} = \nabla \phi

where \\(\\rho\\) is the density, \\(\\mu\\) is the adiabatic compression modulus, \\(\\delta(\\vec{r}-\\vec{r}_s)\\) represents point source, and \\(\\phi\\) is the pressure.

By combining two equations we have 2nd order form:


.. math::

        \rho\frac{\partial^2 \phi}{\partial t^2} - \div\mu\grad\phi = \frac{\partial^2 s}{\partial t^2}\delta(\vec{r}-\vec{r}_s)

By applying Fourier transform to frequency domain we have:

.. math::

        -\rho\omega^2\phi - \div\mu\grad\phi = -\omega ^2 s \delta(\vec{r}-\vec{r}_s)

where \\(\\omega = 2\\pi f\\) is angular frequency. Assuming contant \\(\\mu\\) and \\(\\rho\\) we have Helmholtz equation:

.. math::

        (\grad^2 + k^2)\phi = \mu^{-1}\omega^2 s \delta(\vec{r}-\vec{r}_s)

where \\(\ k = \\omega\\sqrt{\\frac{\\rho}{\\mu}}=\\frac{\\omega}{v} \\) is the wave propagation constant.


Discretization of problem
-------------------------

Artificial Boundary condition
-----------------------------

Sponge boundary
===============

PML boundary
============


.. raw:: html
    :file: examples\refraction.html


Backgrounds
===========


Notebooks
=========

