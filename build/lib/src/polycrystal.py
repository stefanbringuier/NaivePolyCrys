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
        
    def grain_boundary_atoms(self, criteria=0.5):
        """Return indices within ``criteria`` Angstrom of a grain interface.

        Grain tags index the nonperiodic Voronoi centers used by Grains.
        Distances to their separating planes exclude interior atoms without
        comparing atom pairs. Periodic box faces are also candidates because
        pruning uses minimum-image distances. This does not make grain
        generation periodic. Coincident centers and edited geometries that no
        longer follow their grain labels conservatively retain all candidates.
        """
        criteria = float(criteria)
        if not np.isfinite(criteria) or criteria < 0:
            raise ValueError('criteria must be a finite nonnegative distance')
        atoms = self.polycrystal
        if len(atoms) == 0 or criteria == 0:
            return np.empty(0, dtype=int)
        centers = np.asarray(self.crystallitecenters, dtype=float)
        if centers.shape != (self.numcrystallites, 3) or not np.all(np.isfinite(centers)):
            raise ValueError('grain centers must be finite three-dimensional positions')
        tags = atoms.get_tags()
        if np.any(tags < 0) or np.any(tags >= self.numcrystallites):
            raise ValueError('atom tags must index the grain centers')
        positions = atoms.positions
        if not np.all(np.isfinite(positions)):
            raise ValueError('atom positions must be finite')
        lengths = None
        if np.any(atoms.pbc):
            cell = np.asarray(atoms.cell)
            lengths = np.diag(cell)
            if (not np.array_equal(cell, np.diag(lengths))
                    or not np.all(np.isfinite(lengths)) or np.any(lengths <= 0)):
                raise ValueError('periodic pruning requires a positive orthorhombic cell')
        if len(np.unique(tags)) < 2:
            return np.empty(0, dtype=int)
        if lengths is not None:
            periodic = np.flatnonzero(atoms.pbc)
            if np.any(positions[:, periodic] < 0) or np.any(positions[:, periodic] > lengths[periodic]):
                # Wrapped images can overlap far from the unwrapped Voronoi planes.
                return np.arange(len(atoms))

        boundary = np.zeros(len(atoms), dtype=bool)
        for grain in np.unique(tags):
            indices = np.flatnonzero(tags == grain)
            relative = positions[indices] - centers[grain]
            for other in range(self.numcrystallites):
                if other == grain:
                    continue
                delta = centers[other] - centers[grain]
                separation = np.linalg.norm(delta)
                if separation == 0:
                    boundary[indices] = True
                    continue
                distance = separation / 2 - np.dot(relative, delta / separation)
                if np.any(distance < 0):
                    # A moved atom or edited tag invalidates the Voronoi filter.
                    return np.arange(len(atoms))
                boundary[indices] |= distance <= criteria

        if np.any(atoms.pbc):
            for axis in np.flatnonzero(atoms.pbc):
                coordinate = np.mod(positions[:, axis], lengths[axis])
                boundary |= np.minimum(coordinate, lengths[axis] - coordinate) <= criteria
        return np.flatnonzero(boundary)

    def pruneoverlap(self,criteria=0.5,verbose=False):
        """Prune cross-grain pairs strictly closer than ``criteria`` Angstrom.

        Only boundary candidates enter the spatial neighbor search. Process
        original-index pairs in order, keeping the lower index when both atoms
        remain. Atoms already removed cannot cause further deletions. Same-grain
        atoms are preserved; this is geometric cleanup, not physical relaxation.
        """
        boundary = self.grain_boundary_atoms(criteria)
        criteria = float(criteria)
        print("NOTICE: The pruning routine does not account for stoichiometry restrictions.")
        removed = set()
        if len(boundary) > 1:
            candidates = self.polycrystal[boundary]
            first, second = neighborlist.neighbor_list(
                'ij', candidates, criteria, self_interaction=False)
            tags = self.polycrystal.get_tags()
            pairs = sorted({tuple(sorted((int(boundary[i]), int(boundary[j]))))
                            for i, j in zip(first, second)
                            if tags[boundary[i]] != tags[boundary[j]]})
            for i, j in pairs:
                if i not in removed and j not in removed:
                    removed.add(j)
                    if verbose:
                        print("Atom ID %i removed due to overlap!" % j)
        del self.polycrystal[sorted(removed)]
        self.natoms = self.polycrystal.get_global_number_of_atoms()
        self.chemmap = {symbol: i + 1 for i, symbol in
                        enumerate(sorted(set(self.polycrystal.get_chemical_symbols())))}
        self.nspecies = len(self.chemmap)

