from random import randint
from math import acos,pi
import numpy as np
from ase import Atoms

from simbox import OrthoBox
from crystal import Crystal
from polycrystal import Polycrystal

class Grains(OrthoBox):
    def __init__(self,boxsize,ngrains,grainsize=0.5):
        super().__init__(boxsize)
        self.boxsize = boxsize
        self.num = ngrains
        self.size = grainsize
        self.ids = []
        self.graincrystal = {}
        self.natoms = 0
        self.graindist()
        self.graincenters()
        

    def graindist(self,distro='random',dim=3):
        ''' Distrubiton for grain centers '''
        
        if distro == 'random':
            self.points = np.random.random_sample((self.num,dim))
        elif distro == 'normal':
            self.points = np.random.normal(self.size,0.1,(self.num,dim))
             
    def graincenters(self):
        ''' Place grain centers in simulation box '''
        for k in range(3):
            self.points[:,k] *= self.sidelens[k]
        self.points

    def atom_in_grain(self,atompos,grain):
        ''' if atom is closest to this grain center keep'''
        flag = True
        dist_square = np.sum(np.square(atompos-grain[1]))
        for gid,gpos in enumerate(self.points):
            if gid == grain[0]: continue
            tmp_dist_square = np.sum(np.square(atompos-gpos))
            if tmp_dist_square < dist_square:
                flag = False
        return flag
        
    def makepolycrystal(self,names,prunecriteria=1.0):
        ''' Method for creating grains/crystallites in a simulation box resulting in a
        polycrystal configuration.
        
        Input:
           names - list of material/phase for each grain (see src/materials.py)
        '''
        numgrains = self.num

        assert numgrains == len(names)
        graintypes = {gid:names[gid] for gid in range(numgrains)}
        
        polycrystal = Polycrystal(self.boxsize,numgrains,names,self.points)
        polycrystal.setgraintypes(graintypes)
        
        for gid,gcenter in enumerate(self.points):
            print("-------------------------")
            print("Grain id: ", gid)
            self.ids.append(gid)

            # Naive approach: Generate mask crystal such that its bigger than the simulation
            # box. Then calculate mask crystal center of mass and shift mask so its center
            # is at the center of the grain
            
            mask = Crystal(names[gid])
            self.graincrystal[gid] = mask
            ax,ay,az = mask.crystal.cell
            nx,ny,nz = 3*round(self.lx/ax[0]),3*round(self.ly/ay[1]),3*round(self.lz/az[2])
            mask.expandcrystal(nx=nx,ny=ny,nz=nz)

            chemsymbols = mask.crystal.get_chemical_symbols()
            mapchemsymtoid = {k:i+1 for i,k in enumerate(set(chemsymbols))}
            print("-------chemical:id--------")
            for s in mapchemsymtoid.keys():
                print("%s : %i" %(s,mapchemsymtoid[s]))
            print("-------------------------")

            comdiff = gcenter-mask.crystal.get_center_of_mass()
            mask.crystal.positions[:,0] += comdiff[0]
            mask.crystal.positions[:,1] += comdiff[1]
            mask.crystal.positions[:,2] += comdiff[2]
            seed = randint(1,9)*0.1
            mask.xrotate(seed=seed)

            keep = [] #Atom index 
            for atom in range(mask.natoms):
                xyz = mask.crystal.positions[atom,:]
                #check if atom is in box
                if xyz[0] < 0.0 or xyz[0] > self.lx:
                    continue
                elif xyz[1] < 0.0 or xyz[1] > self.ly:
                    continue
                elif xyz[2] < 0.0 or xyz[2] > self.lz:
                    continue

                if self.atom_in_grain(xyz,(gid,gcenter)) == True:
                    keep.append(atom)

            print("Number of atoms in cell: ", len(keep))
            print("Center of grain: ", self.points[gid])
            print("-------------------------")

            polycrystal.embedatoms(mask.crystal,gid,keep)
            
        polycrystal.compress()
        polycrystal.pruneoverlap(criteria=prunecriteria)
        return polycrystal
