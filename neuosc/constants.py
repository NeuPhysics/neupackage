"""
======================================================
Neutrino Physics Constants (:mod:`neuosc.constants`)
======================================================

.. currentmodule:: neuosc.constants

This is a collection of physical constants related to neutrino physics. The overal framework of this module is grabbed from scipy.constants.



Neutrinos Constants
----------------------

============================  =================================================================
``hbarc``                     :math:`\hbar\cdot c = 197.3`
``mass_typical``              the typical mass from experiments = 2eV
``sin_squared_theta_12``      :math:`\sin^2 (\\theta_{12})=0.304`
``delta_mass_squared_21``     :math:`\delta m_{21}^2 = 7.53 \\times 10^{-5} \mathrm{eV}^2`
``sin_squared_theta_23``      :math:`\sin^2 (\\theta_{23})=0.514` for NH
``sin_squared_theta_23_nh``   :math:`\sin^2 (\\theta_{23})=0.514` for NH
``sin_squared_theta_23_ih``   :math:`\sin^2 (\theta_{23})=0.511` for IH
``delta_mass_squared_32``     :math:`\delta m_{32}^2 = 2.44\\times 10^{-3} \mathrm{eV^2}`
``delta_mass_squared_32_nh``  :math:`\delta m_{32}^2=2.44\\times 10^{-3}\mathrm{eV^2}`
``delta_mass_squared_32_ih``  :math:`\delta m_{32}^2 = 2.49\\times 10^{-3}\mathrm{eV^2}`
``sin_squared_theta_13``      :math:`\sin^2 (\\theta_{13})= 2.19\\times 10^{-2}`
============================  =================================================================


:mod:`neuosc.constants` use data from Particle Data Group
2015 Review of Particle Physics Summary Tables [PDGSUMMARY2015]_ .


.. [PDGSUMMARY2015] K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014) and 2015 update.

   http://pdg.lbl.gov/2015/tables/contents_tables.html


Functions
--------------------



"""

######## Physics Constants ########

hbarc = 197.3 # hbar*c in natural units is 1,
# meanwhile, hbar*c = 197.32705 fm*MeV
# which connects length and energy


######## Neutrinos ########

mass_typical = 2 # <2eV


######## Neutrino Mixing ########

sin_squared_theta_12 = 0.304 # $\sin^2 (\theta_{12}) = 0.304 \pm 0.014$

delta_mass_squared_21 = 7.53e-5 # eV^2; $\delta m_{21}^2$

sin_squared_theta_23 = sin_squared_theta_23_nh = 0.514 # for normal hierarchy;
# $\sin^2 (\theta_{23})$
# the naming is meant to use normal hierarchy as the default for 23
# however this is just a convention here.
# ------------------------------------------------

sin_squared_theta_23_ih = 0.511 # for inverted hierarchy;
# $\sin^2 (\theta_{23})$
# ------------------------------------------------

delta_mass_squared_32 = delta_mass_squared_32_nh = 2.44e-3 # eV^2; $\delta m_{32}^2$
delta_mass_squared_32_ih = 2.49e-3 # eV^2; $\delta m_{32}^2$

sin_squared_theta_13 = 2.19e-2 # $\sin^2 (\theta_{13})$




# Conversion Between SI and Natural Unit


def mev2km(x):
    """
    Find the corresponding km length scale of x MeV

    **Parameters**


    x : real number
        How many MeVs

    **Returns**


    return : return the converted value

    **Examples**

    >>> mev2km(1)


    """
    # Convert

    return 1
