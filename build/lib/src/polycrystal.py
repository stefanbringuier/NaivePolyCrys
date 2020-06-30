#Dependencies
import numpy as np
from ase import Atoms
from ase import neighborlist

#Internal
from simbox import OrthoBox

class Polycrystal(OrthoBox):
    def __init__(self,boxsize,numcrystallites,crystallitenames,crystallitecenters):
        super().__init__(boxsize)
        self.numcrystallites = numcrystallites
        self.crystaltype = crystallitenames
        self.crystallitecenters = crystallitecenters
        self.crystallites = [Atoms() for n in range(numcrystallites)]
        self.crystalliteids = [[i] for i in range(numcrystallites)]

    def setgraintypes(self,typemap):
        self.graintypes = typemap
        
    def embedatoms(self,maskcrystal,grainid,keepatoms,compress=True):
        crystal = maskcrystal.copy()
        slicecrystal = crystal[keepatoms]
        natoms = slicecrystal.get_global_number_of_atoms()
        self.crystalliteids[grainid] = [grainid for i in range(natoms)]
        self.crystallites[grainid] =  slicecrystal
            
    def compress(self):
        polycrystalbox = self.sidelens+self.angles
        self.polycrystal = Atoms(cell=polycrystalbox,pbc=True)
        for i,c in enumerate(self.crystallites):
            cids = self.crystalliteids[i]
            c.set_tags(cids)
            self.polycrystal += c
        self.natoms = self.polycrystal.get_global_number_of_atoms()
        chemsym = set(self.polycrystal.get_chemical_symbols())
        self.chemmap = {c:i+1 for i,c in enumerate(chemsym)}
        self.nspecies = len(chemsym)
        
    def pruneoverlap(self,criteria=0.5,verbose=False):
        ''' TODO: remove atoms overlapping assume PBC conditions'''
        print("NOTICE: The pruning routine does not accoutn for stoichiometry restrictions.")
        
        cutoff = neighborlist.natural_cutoffs(self.polycrystal)
        neighbors = neighborlist.NeighborList(cutoff,
                                                 self_interaction=False,
                                                 bothways=False)
        neighbors.update(self.polycrystal)

        isremoved = []
        for a in self.polycrystal:
            ineigh = neighbors.get_neighbors(a.index)
            for j in ineigh[0]:
                r = self.polycrystal.get_distance(a.index,j,mic=True)
                if r < criteria and j not in isremoved:    
                    isremoved.append(j)
                    if verbose:
                        print("Atom ID %i removed due to overlap!" %(j))

        del self.polycrystal[isremoved]
        self.natoms = self.polycrystal.get_global_number_of_atoms()
        #NOTE! Assuming the different species types remains the same
