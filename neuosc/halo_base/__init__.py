from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

omegav = 1

mu = 50 * omegav # in unit of omegav where omegav = 1
refl = 0.1 # reflected neutrinos
thetav = np.arcsin( np.sqrt( 0.3 ) )


## Prepare

# Initial set up
rhoF0 = np.array( [ [ 1, 0 ], [ 0, -1 ] ] )/2
rhobB0 = np.array( [ [ 0, 0 ], [ 0, 0 ] ] )/2

# Vacuum Hamiltonian in unit of omega
hvac = np.array( [ [ -np.cos(thetav), np.sin(thetav) ], [ np.sin(thetav), np.cos(thetav) ] ] )



def hamilF( rhoF, rhobB, muValue ):
    """
    Calculate the hamiltonian for forward beams
    """

    hnunuF = 2 * muValue * ( rhoF - rhobB )

    return hvac + hnunuF


def hamilbB( rhobB, rhoF, muValue ):

    return -hvac + 2 * muValue * ( rhoF - rhobB )



## Define Vacuum Solver

def vacSolver(N, dx, initx, initrho, hamil = None):

    dx = np.abs(dx)

    if hamil == None:
        hamil = hvac

    rhoa = np.zeros( N )
    rhobr = np.zeros( N )
    rhobi = np.zeros( N )

    rhoa[0], rhobr[0], rhobi[0] = initrho[0,0].real, initrho[0,1].real, initrho[0,1].imag


    for i in np.arange(0,N-1):

        fdarhs = - hamil[0,1].real * rhobi + hamil[0,1].imag * rhobr
        fdbrrhs = - hamil[0,1].imag * rhoa + hamil[0,0].real * rhobi
        fdbirhs = - hamil[0,0].real * rhobr + hamil[0,1].real * rhoa


        rhoa[i+1] = rhoa[i] +  dx * fdarhs[i]
        rhobr[i+1] = rhobr[i] +  dx * fdbrrhs[i]
        rhobi[i+1] = rhobi[i] +  dx * fdbirhs[i]

    rhoarray = np.array([ [ [ rhoa[i], rhobr[i] + 1j *  rhobi[i] ], [ rhobr[i] - 1j *  rhobi[i], -rhoa[i] ]  ] for i in np.arange(0,N-1)])

    return rhoarray



# Solve Vacuum Problem


def initL(L, dz):

    Nz = int(L/dz)

    rhoFinitL = vacSolver(Nz, dz, 0, rhoF0 )

    # determine the flavor at z = L

    rhoL = rhoFinitL[-1]
#     print(len(rhoFinitL))

    probFinitL = [2*rhoFinitL[i,0,0].real for i in np.arange( len(rhoFinitL) )]

    return {'rhoL': rhoL, 'prob': probFinitL, 'rho': rhoFinitL  }




## Beam Solver

def BeamSolver(N, dx, initx, initrho, hamil):

    dx = np.abs(dx)

    rhoa = np.zeros( N )
    rhobr = np.zeros( N )
    rhobi = np.zeros( N )

    rhoa[0], rhobr[0], rhobi[0] = initrho[0,0].real, initrho[0,1].real, initrho[0,1].imag

    for i in np.arange(0,N-1):

        fdarhs = - hamil[i,0,1].real * rhobi + hamil[i,0,1].imag * rhobr
        fdbrrhs = - hamil[i,0,1].imag * rhoa + hamil[i,0,0].real * rhobi
        fdbirhs = - hamil[i,0,0].real * rhobr + hamil[i,0,1].real * rhoa


        rhoa[i+1] = rhoa[i] +  dx * fdarhs[i]
        rhobr[i+1] = rhobr[i] +  dx * fdbrrhs[i]
        rhobi[i+1] = rhobi[i] +  dx * fdbirhs[i]

    rhoarray = np.array([ [ [ rhoa[i], rhobr[i] + 1j *  rhobi[i] ], [ rhobr[i] - 1j *  rhobi[i], -rhoa[i] ]  ] for i in np.arange(0,N-1)])

    return rhoarray


# Calculate hamiltonian

def BeamHamil(rhoF, rhoB, mu):

    return 2 * mu * ( rhoF - rhoB )
