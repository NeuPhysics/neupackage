# This is a script that solves the linear stability for MZA solutions for box spectrum.

# Add current direction to sys.path
import os, sys
cbd = os.getcwd()
sys.path.insert(0, cbd+'/python/dispersion-relation/')

from __future__ import division
import numpy as np
import sympy as sp
import mpmath as mp

mp.dps = 15
mp.pretty = True

# import the spectra as parameters
import boxspectra


# Define function of MZA solutions


def IntFun0n(n, ct1, ct2):
    return (-np.log( (1- ct1 * n)/(1-ct2*n) ) ) / n

def IntFun1n(n, ct1, ct2):
    return ((-ct1 + ct2) +np.log( (ct2*n-1)/(ct1*n-1) ) / n ) / n

def IntFun2n(n, ct1, ct2):
    return (-(ct1-ct2) * ( (ct1+ct2) +2/n )+ 2* np.log( (ct2 * n-1)/(ct1 * n-1) ) /n**2)/(2 * n)



for i in range( len(spectC1) ):
    print spectC1[i][0][0], spectC1[i][0][1]

IntFun0n(1, spectC1[0][0][0], spectC1[0][0][1])

def ConAxialSymOmegaNMZApEqnLHSComplex(omega,k, spect):

    spectM = spect;

    IntFun0nFFM = Total[#[[2]] IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
    IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
    IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

    eqnLHSM = omega-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
    return eqnLHSM


ConAxialSymOmegaNMZApEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;


IntFun0nFFM = Total[#[[2]] IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omega-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);





eqnLHSM



]
