import numpy as np
from ase.spacegroup import crystal
from ase.build import make_supercell
from materials import *

class Crystal:
    """ Interface with ASE crystal builder to generate crystal. Constructs
    an crystal object.

    Input: crystalname -> material phase.

    """
    def __init__(self,crystalname):
        self.name = crystalname
        self.units = {'lengths':'Angstroms',
                      'angles': 'Degrees'}
        self.length = []
        self.angles = []
        self.cell = []
        self.species = []
        self.basis = []
        self.spacegroup = None
        self.getcrystal()
        
    def getcrystal(self):
        self.getsetting()
        self.crystal = crystal(self.species,
                       basis=self.basis,
                       spacegroup=self.spacegroup,
                       cellpar = self.cell)
        self.natoms = self.crystal.get_global_number_of_atoms()
        
    def expandcrystal(self,nx=1,ny=1,nz=1):
        P = [[nx,0,0],
             [0,ny,0],
             [0,0,nz]]
        self.crystal = make_supercell(self.crystal,P)
        self.natoms = self.crystal.get_global_number_of_atoms()
        
    def getsetting(self):
        setting = materials[self.name]
        self.ref = setting['ref']
        self.lengths = setting['lengths']
        self.angles = setting['angles']
        self.cell = self.lengths + self.angles
        self.species = setting['species']
        self.basis = setting['basis']
        self.spacegroup = setting['spacegroup']
        

    def xrotate(self,seed=1):
        """
        Use Rodrigues rotation formula axis x vector = matrix
        then use Rodrigues to rotate vector about axis for theta

        Input: 
        crystal - perfect crystal
        seed - seed for angle (unitless)
        Output:
        rotcrys - rotated crystal with seed about random axis
        """
        atompos = self.crystal.positions
    
        angle = seed*np.pi/2. #radians
        axis = np.random.random_sample((3)) #Unit vector


        axisXvec = np.array([[0.0, -axis[2], axis[1]],
                         [axis[2], 0.0 , -axis[0]],
                         [-axis[1], axis[0], 0.0]])

        #Rodrigues Formula
        rotmat = (np.identity(3) + axisXvec * np.sin(angle) 
              + np.dot(axisXvec,axisXvec) * (1. - np.cos(angle)))


        #Rotate Crystal
        natoms = self.crystal.get_global_number_of_atoms()
        rotcrys = np.ndarray((natoms,3),dtype=float)

        for i in range(natoms):
            rotcrys[i,:] = np.dot(rotmat,atompos[i,:])

        self.crystal.set_positions(rotcrys)
