# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 20:36:13 2023

@author: Ajay Khanna
@Lab: Dr. Isborn's Lab
@Place: UC Merced
"""
import numpy as np
import random
import matplotlib.pyplot as plt

def pairwise_potential(epsilon, sigma, r_i, r_j, r_c):
    """
    It calculates the energy of interaction between two atoms, given their positions and the cut-off
    radius
    
    :param epsilon: Control the strenght of the interaction
    :param sigma: Length, different for different atoms
    :param r_i: Coordinates of atom i
    :param r_j: Coordinates of atom j
    :param r_c: cut-off radius, beyond this interaction dies to zero
    :param r_c = 2^(1/6) * sigma # For Argon
    :return: The energy of the interaction between two atoms
    """

    if abs(r_j - r_i) >= r_c:
        energy = 0
    else:
        repulsion = np.power((sigma/abs(r_i - r_j)), 12)
        attaction = np.power((sigma/abs(r_i - r_j)), 6)
        energy = 4 * epsilon * (repulsion - attaction) + epsilon
    
    return energy

if __name__ == '__main__':
    natoms = int(10)     # Total Argon Atoms
    sigma = float(3.405) # Argon
    r_c = float(2.0**(1/6) * sigma) # Original paper
    epsilon = float(0.238) # Argon
    coordinates = random.sample(range(-10,10), natoms) # Random Coordinates
    lj_interactions = []
    for i in range(natoms):
        for j in range(natoms):
            if i<j: # Make sure no self-interaction or double counting
                lj_interactions.append(pairwise_potential(epsilon, sigma,
                                                    coordinates[i], coordinates[j],
                                                    r_c))
    
    print(f'Pairwise LJ Potential Energy: {lj_interactions}')
    plt.xlabel('#Pairwise Potential')
    plt.ylabel('LJ Potential Energy (arb. units)')
    plt.plot(np.arange(len(lj_interactions)), lj_interactions)
