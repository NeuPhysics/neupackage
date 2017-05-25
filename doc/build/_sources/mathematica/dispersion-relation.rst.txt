Dispersion Relation
============================


The package is located at `mma/dispersion-relation.wl`. The corresponding test notebook is located at `mma/package-test/dispersion-relation-package-test.nb`.

For a reference of paper about this topic, please read

1. Izaguirre, I., Raffelt, G., & Tamborra, I. (2017). `Fast Pairwise Conversion of Supernova Neutrinos: A Dispersion Relation Approach <https://doi.org/10.1103/PhysRevLett.118.021101>`_. Physical Review Letters, 118(2), 21101.

Two Beams
----------------------



N Beams
----------------------


Box Spectra
-------------------------

In the package, a box spectrum is defined as

.. math::
   \{ \{ \{ u_1, u_1' \}, g_1  \} ,\{ \{ u_2, u_2' \}, g_2  \}, ,\{ \{ u_3, u_3' \}, g_3  \} , \cdots \},

where :math:`u_i` is the start :math:`\cos\theta_i` value and :math:`u_i'` is the ending value of :math:`\cos\theta_i`, within these two values, we have the spectrum value :math:`g_1`.

The functions defined in this section can take in spectrum of arbitrary segments.

During the calculation of any quantities in this problem, the integral

.. math::
   I_m = \int_{c_1}^{c_2} \frac{u^m}{1-n u} \, du

is widely used. These integrals can be calculated analytically.

.. function:: IntFun0n[n,c1,c2]

   Calculates the value of :math:`I_0` for given :math:`n`, :math:`c_1`, and :math:`c_2`.

   :param n: the variable :math:`n`
   :param c1: the lower limit of the integral
   :param c2: the upper limit of the integral
   :rtype: a real or complex number



.. function:: IntFun1n[n,c1,c2]

   Calculates the value of :math:`I_1` for given :math:`n`, :math:`c_1`, and :math:`c_2`.

   :param n: the variable :math:`n`
   :param c1: the lower limit of the integral
   :param c2: the upper limit of the integral
   :rtype: a real or complex number



.. function:: IntFun2n[n,c1,c2]

   Calculates the value of :math:`I_2` for given :math:`n`, :math:`c_1`, and :math:`c_2`.

   :param n: the variable :math:`n`
   :param c1: the lower limit of the integral
   :param c2: the upper limit of the integral
   :rtype: a real or complex number


.. function:: ConAxialSymOmegaNMAA[n,spect_optional]

   Calculates :math:`\omega(n)` for MAA solution for given spectrum.

   :param n: the variable :math:`n`
   :param spect: the input spectrum, which is optional. The default spectrum is {{{0.3,0.6},1},{{0.6,0.9},-1}}
   :rtype: real or complex number

.. function:: ConAxialSymOmegaNMZA[n,spect_optional]

   Calculates :math:`\omega(n)` for MZA solutions for given spectrum.

   :param n: the variable :math:`n`
   :param spect: the input spectrum, which is optional. The default spectrum is {{{0.3,0.6},1},{{0.6,0.9},-1}}
   :rtype: a list of real or complex number, {MZA+, MZA-}
