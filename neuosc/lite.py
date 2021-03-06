"""
======================================================
NeuOsc Lite (:mod:`neuosc.lite`)
======================================================

.. currentmodule:: neuosc.lite


A very lite version of this package with all frequently
used functions and files in one place.


Classes and Functions
------------------------



"""

# A module for neutrino oscillations

import numpy as np

## Define neutrino mixing matrix which brings the vacuum energy basis to flavor basis.


# class par(object):





def pauli_matrices(n):
    """
    Pauli Matrices

    **Parameters**

    n : 0,1,2,3
        specify which matrix to use.
        n = 0: identity
        n = 1: sigma_1
        n = 2: sigma_2
        n = 3: sigma_3

    **Returns**

    returns the corresponding Pauli matrix as a np.array()

    **Notes**

    n = 0 is identity.

    **Usage**

    - pauli_matrices(n)  is :math:`\sigma_n`;
    - n runs from 0 to 3, where pauli_matrices(0) is identity matrix.

    **Examples**

    >>> from neuosc.lite import *
    >>> pauli_matrices(1)
    >>> [[0,1.],[1,0.]]

    """

    matrices = np.array([ [[1.0,0], [0,1.0]], [[0,1.],[1,0.]], [[0,-1.j],[1.j,0]], [[1.,0],[0,-1.]] ])

    return matrices[n-1]




## Mixing Matrix for Two Flavor Oscillation

def mixing2():  ### 2 flavor case with fixing mixing angles from recent results
    """neutrino mixing matrix with : from vacuum basis to flavor basis"""

    sin2thetav = 0.917
    cos2thetav = 0.4
    return np.array([[cos2thetav,sin2thetav],[-sin2thetav,cos2thetav]])

def mix2(thetav):   ### 2 flavor case
    """neutrino mixing matrix: from """

    return np.array([[np.cos(2*thetav),np.sin(2*thetav)],[-np.sin(2*thetav),np.cos(2*thetav)]])




class Neutrino(object):
    """
    Neutrino class which contains energy, mixing angle, mass differences and vacuum oscillation frequency.
    """

    def __init__(self, energy, mixing_angle, mass_difference):
        """Initialize a new neutrino object"""
        self.energy = energy
        self.mixing_angle = mixing_angle
        self.mass_difference = mass_difference

    def mixing_matrix(self):
        """Unitary matrix that rotates the state from vacuum basis to flavor basis"""
        return np.array([ [ np.cos(self.mixing_angle), np.sin(self.mixing_angle) ],[ -np.sin(self.mixing_angle),np.cos(self.mixing_angle) ] ])

    def omega(self):
        """Frequency omega of the neutrino oscillation in vacuum,

        :math:`\omega = \Delta m^2 / 2 E`
        """
        return self.mass_differences/(2*self.energy)



class MSW(object):
    """ MSW effect calculations

    Attributes:
    energy: energy of the neutrino"""

    def __init__(self, energy):
        self.energy = energy

    def adiabatic_p(self, mixing_angle, lambdahat):
        """ Survival probability of adiabatic evolution

        Attributes:
        
        - mixing_angle: The mixing angle for two flavor oscillation;
        - lambdahat: :math:`\hat\lambda = \lambda/\omega` where :math:`\omega` is the vacuum frequency."""

        theta = mixing_angle
        cos2thetam =  ( np.cos(2*theta) - lambdahat ) /( np.sqrt( lambdahat**2 + 1 - 2*lambdahat*np.cos(2*theta) ) )

        return (1 + np.cos(2*theta)*cos2thetam )/2

    def solar_landau_zener_p(self, mixing_angle, lambdahat):
        """ Survival probability with Landau Zener transition

        Attributes:

        - mixing_angle: The mixing angle for two flavor oscillation;
        - lambdahat: :math:`\hat\lambda = \lambda/\omega` where :math:`\omega` is the vacuum frequency."""

        theta = mixing_angle
        cos2thetam =  ( np.cos(2*theta) - lambdahat  )/( np.sqrt( lambdahat**2 + 1 - 2*lambdahat*np.cos(2*theta) ) )

        gamma = (2.554*10**(3))*( (np.sin(2*theta) )**2 )/( lambdahat*np.cos(2*theta) )
        pf = np.exp(-np.pi*gamma/2)

        return (1 + (1 - 2*pf)*np.cos(2*theta)*cos2thetam )/2








# End of module
