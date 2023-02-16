import numpy as np
from numpy import newaxis
import matplotlib.pyplot as plt

# Physical Constant
avogadro_constant = 6.02214086e23
boltmann_constant = 1.38064852e-23

def check_wall(pos, vels, box):
    """
    > If a particle hits the wall, its velocity is updated in the opposite direction.
    > Enforcing reflective bounday condistions
    
    The function takes in three arguments: 
    1. `pos`: the position of the particles
    2. `vels`: the velocity of the particles
    3. `box`: the box dimensions
    
    The function is written in a way that it can be used for any number of particles and returns nothing 
    :param pos: the position of the particles
    :param vels: the velocities of the particles
    :param box: the box size, in this case, it's a 2D box with size (0,1) in the x-direction and (0,2)
    in the y-direction
    """
    ndims = len(box)
    for i in range(ndims):
        vels[((pos[:,i] <= box[i][0]) | (pos[:,i] >= box[i][1])), i] *= -1

def intergrate(pos, vels, forces, mass, dt):
    """
    Simple forward Euler integrator, propogating system in time

    It takes the current position, velocity, force, and mass of each particle, and uses them to
    calculate the new position and velocity of each particle
    
    :param pos: the positions of the particles
    :param vels: the velocities of the particles
    :param forces: a 2D array of shape (n_particles, n_dimensions)
    :param mass: mass of each particle
    :param dt: time step
    """
    
    pos += vels * dt
    vels += (forces * dt) / mass[np.newaxis].T
        
def computeForces(mass, vels, temp, relax, dt):
    """
    It computes the Langevin force on all atoms
    
    :param mass: The mass of each atom
    :param vels: the velocities of the atoms
    :param temp: The temperature of the system
    :param relax: relaxation time
    :param dt: The time step
    :return: The force on each atom in each dimension.
    """    
    natoms, ndims = vels.shape
    sigma = np.sqrt((2.0 * mass * temp * boltmann_constant) / (relax * dt))
    noise = np.random.randn(natoms, ndims) * sigma[np.newaxis].T
    force = - ((vels * mass[newaxis].T) / relax) + noise
    
    return force


def run_reflectiveBC_MD(**args):
    
    natoms, box, dt, temp = args['natoms'], args['box_dim'], args['dt'], args['temp']
    mass, relax, nsteps = args['mass'], args['relax'], args['nsteps']
    ofname, freq, radius = args['ofname'], args['freq'], args['radius']

    dim = len(box)
    pos = np.random.rand(natoms, dim)
    
    for i in range(dim):
        pos[:,i] = box[i][0] + (box[i][1] - box[i][0]) * pos[:,i]
    
    vels = np.random.rand(natoms, dim)
    mass = np.ones(natoms) * mass / avogadro_constant
    radius = np.ones(natoms) * radius
    step = 0
    
    output = []
    
    while step <= nsteps:
        step +=1
        
        #Compute all forces
        forces = computeForces(mass, vels, temp, relax, dt)
        
        # Propagate the system in time
        intergrate(pos, vels, forces, mass, dt)
        
        # Check if any particle has hit the wall
        check_wall(pos, vels, box)
        
        # Log output
        ins_temp = np.sum(np.dot(mass, (vels - vels.mean(axis=0))**2)) / (boltmann_constant * dim * natoms)
        output.append([step * dt, ins_temp])
        
        #if not step % freq:
        #    dump.writeOutput(ofname, natoms, step, box, radius=radius, v=vels)
     
    return np.array(output)
    
if __name__ == '__main__':

    params = {
        'natoms': 1000,
        'temp': 300,
        'mass': 0.001,
        'radius': 120e-12,
        'relax': 1e-13,
        'dt': 1e-15,
        'nsteps': 10000,
        'freq': 100,
        'box_dim': ((0, 1e-8), (0, 1e-8), (0, 1e-8)),
        'ofname': 'traj-hydrogen-3D.dump'
    }

    output = run_reflectiveBC_MD(**params)

    plt.plot(output[:,0] * 1e12, output[:,1])
    plt.xlabel('Time (ps)')
    plt.ylabel('Temp (K)')
    plt.show()
