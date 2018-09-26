# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 12:41:17 2018

@author: Haldrian
"""

import math
import sys

def computeN (N0, f, dl, nGrid):
# Compute the number of metal ions in solution
########################################################################################
    N = -0.5*(f[0] + f[nGrid-1]) #summing up the metal ions in solution
    for i in range(nGrid):
        N += f[i]
    N *= dl*N0 
    return N

def computeV (N0, beta, f, nGrid, dl, v):
# Compute the electrostatic potential in the solution
########################################################################################
    sigma = [0.0 for i in range(nGrid)] #create a list
    for i in reversed(range(nGrid-1)):
        sigma[i] = sigma[i+1] + 0.5*(f[i] + f[i+1])*dl
#
    v[0] = 2*N0*sigma[0] #where from?
    for i in range(1,nGrid):
        v[i] = v[i-1] + N0*beta*0.5*(sigma[i-1] + sigma[i])*dl

def addBias (beta, vG, nGrid, l, v):
# Compute the change in electrostatic potential in the solution due to a bias
########################################################################################
    for i in range(nGrid):
        v[i] += (3.0 + l[i]*beta)*vG/(6.0 + beta)

def dissolve():
    ''' This function solves the equation of motion for the dissolution of an electrode.'''
# Set the parameters that define the simulation.
# nSteps The number of
# timeStep The time step for the update
# nGrid The number of points on the grid
# N0 The number of metal ions corresponding to the reference concentration
# gamma The solvation energy
# theta The temperature
# beta The strength of the electrostatic interaction
########################################################################################
    Ha = 27.211386 #
    a0 = 0.529177210
#
    z = 2.0
    epsilon = 80.0
    b = 5.0
    W = 100.0
    A = 1.0*W*W
    c0 = 6e-4 #mol per litre, atoms per angstrom3
    kT = 0.025
    g0 = -4.8     
    #
    nSteps = 1000000
    timeStep = 1e-5
    nGrid = 1001
    U = 2*3.1416*b*a0*z*z*Ha/(3.0*A)
    N0 = c0*A*W
    gamma = -g0/U
    theta = kT/U
    beta = 6*W/(b*epsilon)
    
    # Set the bias flag. 0 ==> Open circuit; 1 ==> Bias applied
########################################################################################
    biasFlag = 0
    vBias = 0.0
    muRL = 0.0

    # Build the grid
########################################################################################
    dl = 1.0/(nGrid-1)
    l = [i*dl for i in range(nGrid)] #creates coordinates along the solution axes with gaps dl
                                    # l being x axis
                                    
    # Initialise the concentration and the electrostatic potential
########################################################################################
    f = [1.0e-10 for i in range(nGrid)]
    v = [0.0 for i in range(nGrid)]
    s = 0
    # Evolve the concentration
########################################################################################
    
    while True:
        '''line 97-150 deleted because code is still processed for publication'''
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
# Write out the final concentration profile
########################################################################################
    N = computeN(N0, f, dl, nGrid)
    f0 = [math.exp(-(v[i] - gamma)/theta) for i in range(nGrid)]
    filename = 'steps' + str(nSteps) + ' kT = ' + str(kT)
    outFile = open('%s.csv' % filename, 'w')
    outFile.write('N, U (eV), N0, gamma, theta, beta\n')
    outFile.write(str(N) + ', ' + str(U) + ', ' + str(N0) + ', ' + str(gamma) + ', ' + str(theta) + ', ' + str(beta) + '\n')
    outFile.write('x, f, f0, v-v0\n')
    for i in range(nGrid):
        outFile.write(str(l[i]) + ', ' + str(f[i]) + ', ' + str(f0[i]) + ', ' + str(v[i]-v[0]) + '\n')
    outFile.close()
    print("FOUND IT \n")
    print('steps taken =' + str(s))
    
# Execute the main code if run as a script.
if __name__ == "__main__":
    dissolve()