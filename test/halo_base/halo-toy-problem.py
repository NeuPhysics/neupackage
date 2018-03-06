from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

output = "prob-test.csv"
# Record probability for every rstep element
rstep = 100

# Number of iterations
Ntop = 100



## Setup Parameters

# delta z
dz = 0.0001

# Length of z direction
L = 1

# x range

zarray = np.arange( 0, L, dz )

# Physical Params
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



def hamilF( rhoB, muValue ):
    """
    Calculate the hamiltonian for forward beams
    """

    hnunuF = 2 * muValue * ( rhoB )

    return hvac + hnunuF


def hamilB( rhoF, muValue ):

    return -hvac + 2 * muValue * ( rhoF )



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


def initL():

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

def BeamHamil(rho, mu):

    return 2 * mu * ( rho )



############################
## Solve the problem
#############################

N = int( L/dz )

# rhoFIter = np.zeros( (2, N) )
# rhoBIter = np.zeros( (2, N) )

rhoFIter = initL()['rho']
rhoBIter = vacSolver( N, -dz, L, initL()['rhoL'] )
probStored = [ [rhoFIter[i,0,0].real for i in np.arange(0,N-1,rstep)], [rhoBIter[i,0,0].real for i in np.arange(0,N-1,rstep)] ]

# print( np.array(probStored).shape )

for i in tqdm(np.arange(Ntop)):

    hamilBIter = BeamHamil(rhoFIter[::-1], mu )
    hamilFIter = BeamHamil(rhoBIter[::-1], mu*refl )


    rhoFIter = BeamSolver( int(L/dz), dz, 0, rhoF0, hamilFIter)

    rhoB0 = rhoFIter[-1]

    rhoBIter = BeamSolver( int(L/dz), -dz, L, rhoB0, hamilBIter)

    probFIter = np.array([rhoFIter[i,0,0].real for i in np.arange(0,N-1,rstep)])
    probBIter = np.array([rhoBIter[i,0,0].real for i in np.arange(0,N-1,rstep)])

    if (i % ( int(Ntop/10) ) == 0) | (np.abs(Ntop - i)<5 ):

        probStored.append( probFIter )
        probStored.append( probBIter )
        # print(i)
        # plt.plot( probFIter, label='Forward' )
        # plt.plot( probBIter[::-1], label='Backward' )
        # plt.legend()
        # plt.show()

# print( np.array(probStored).shape )
np.savetxt(output, probStored, delimiter=",")
