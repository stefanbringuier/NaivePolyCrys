"""Geometric boundary and overlap checks; no force fields or dynamics."""
import contextlib
import io
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from polycrystal import Polycrystal


def sample(positions, tags, centers=((2., 5., 5.), (8., 5., 5.)),
           symbols=None, pbc=False):
    value = Polycrystal(np.diag([10., 10., 10.]), len(centers),
                        ['Si'] * len(centers), np.array(centers))
    value.polycrystal = Atoms(symbols or ['Si'] * len(positions),
                             positions=np.array(positions).reshape(-1, 3),
                             tags=tags, cell=[10., 10., 10.], pbc=pbc)
    value.natoms = len(positions)
    return value


def prune(value, criteria):
    with contextlib.redirect_stdout(io.StringIO()):
        value.pruneoverlap(criteria)


def brute_pairs(value, criteria):
    atoms = value.polycrystal
    tags = atoms.get_tags()
    return [(i, j) for i in range(len(atoms)) for j in range(i + 1, len(atoms))
            if tags[i] != tags[j] and atoms.get_distance(i, j, mic=True) < criteria]


class GrainBoundaryTests(unittest.TestCase):
    def test_tracked_build_copy_matches_source(self):
        root = Path(__file__).resolve().parents[1]
        self.assertEqual((root / 'src/polycrystal.py').read_bytes(),
                         (root / 'build/lib/src/polycrystal.py').read_bytes())

    def test_internal_boundary_excludes_grain_interiors(self):
        value = sample([[2, 5, 5], [4.8, 5, 5], [5.2, 5, 5], [8, 5, 5]],
                       [0, 0, 1, 1])
        np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [1, 2])

    def test_pruning_preserves_close_same_grain_atoms(self):
        value = sample([[4.8, 5, 5], [4.9, 5, 5], [5.1, 5, 5]], [0, 0, 1])
        prune(value, .5)
        np.testing.assert_allclose(value.polycrystal.positions, [[4.8, 5, 5], [4.9, 5, 5]])
        self.assertEqual(value.natoms, 2)

    def test_user_cutoff_not_limited_by_covalent_radii(self):
        value = sample([[4., 5, 5], [6., 5, 5]], [0, 1], symbols=['H', 'H'])
        prune(value, 2.1)
        self.assertEqual(value.natoms, 1)

    def test_periodic_boundary_candidates_and_pruning(self):
        value = sample([[.1, 5, 5], [9.9, 5, 5], [2, 5, 5], [8, 5, 5]],
                       [0, 1, 0, 1], pbc=True)
        np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [0, 1])
        prune(value, .5)
        self.assertEqual(value.natoms, 3)
        self.assertEqual(brute_pairs(value, .5), [])

    def test_nonperiodic_box_faces_are_not_grain_boundaries(self):
        value = sample([[.1, 5, 5], [9.9, 5, 5]], [0, 1])
        np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [])
        prune(value, .5)
        self.assertEqual(value.natoms, 2)

    def test_exact_cutoff_is_not_an_overlap(self):
        value = sample([[4.75, 5, 5], [5.25, 5, 5]], [0, 1])
        prune(value, .5)
        self.assertEqual(value.natoms, 2)

    def test_deleted_atom_does_not_remove_another_atom(self):
        centers = ((5, 5, 5),) * 3
        value = sample([[5, 5, 5], [5.4, 5, 5], [5.8, 5, 5]], [0, 1, 2], centers)
        prune(value, .5)
        np.testing.assert_allclose(value.polycrystal.positions, [[5, 5, 5], [5.8, 5, 5]])

    def test_coincident_centers_are_conservative(self):
        value = sample([[2, 5, 5], [2.1, 5, 5]], [0, 1], ((2, 5, 5),) * 2)
        np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [0, 1])
        prune(value, .5)
        self.assertEqual(value.natoms, 1)

    def test_only_boundary_subset_reaches_neighbor_search(self):
        value = sample([[2, 5, 5], [4.8, 5, 5], [5.2, 5, 5], [8, 5, 5]], [0, 0, 1, 1])
        from polycrystal import neighborlist
        original = neighborlist.neighbor_list
        seen = []
        def recording(quantities, atoms, cutoff, **kwargs):
            seen.append(atoms.positions.copy())
            return original(quantities, atoms, cutoff, **kwargs)
        with patch.object(neighborlist, 'neighbor_list', side_effect=recording):
            prune(value, .5)
        self.assertEqual(len(seen), 1)
        np.testing.assert_allclose(seen[0], [[4.8, 5, 5], [5.2, 5, 5]])

    def test_empty_single_grain_and_zero_cutoff_are_noops(self):
        values = [sample([], []), sample([[2, 5, 5], [2.1, 5, 5]], [0, 0])]
        for value in values:
            original = value.polycrystal.positions.copy()
            np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [])
            prune(value, .5)
            np.testing.assert_array_equal(value.polycrystal.positions, original)
        value = sample([[5, 5, 5], [5, 5, 5]], [0, 1])
        prune(value, 0)
        self.assertEqual(value.natoms, 2)

    def test_invalid_cutoffs_refuse_without_mutation(self):
        for criteria in [-1, float('nan'), float('inf')]:
            value = sample([[4.8, 5, 5], [5.2, 5, 5]], [0, 1])
            original = value.polycrystal.positions.copy()
            with self.assertRaises(ValueError):
                prune(value, criteria)
            np.testing.assert_array_equal(value.polycrystal.positions, original)

    def test_moved_atoms_fall_back_without_losing_an_overlap(self):
        value = sample([[2, 5, 5], [2.1, 5, 5], [8, 5, 5]], [0, 1, 1])
        np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [0, 1, 2])
        prune(value, .5)
        self.assertEqual(value.natoms, 2)
        self.assertEqual(brute_pairs(value, .5), [])

    def test_out_of_box_periodic_images_fall_back(self):
        for second in [12.1, -7.9]:
            centers = ((2, 5, 5), (second - .1, 5, 5))
            value = sample([[2, 5, 5], [second, 5, 5]], [0, 1], centers, pbc=True)
            np.testing.assert_array_equal(value.grain_boundary_atoms(.5), [0, 1])
            prune(value, .5)
            self.assertEqual(value.natoms, 1)

    def test_moved_atom_cannot_bypass_periodic_cell_refusal(self):
        value = sample([[2, 5, 5], [2.1, 5, 5]], [0, 1], pbc=True)
        value.polycrystal.cell[0, 1] = 1
        original = value.polycrystal.positions.copy()
        with self.assertRaisesRegex(ValueError, 'orthorhombic'):
            prune(value, .5)
        np.testing.assert_array_equal(value.polycrystal.positions, original)

    def test_species_metadata_tracks_a_species_removed_entirely(self):
        value = sample([[4.9, 5, 5], [5.1, 5, 5]], [0, 1], symbols=['Si', 'Ge'])
        prune(value, .5)
        self.assertEqual(value.polycrystal.get_chemical_symbols(), ['Si'])
        self.assertEqual(value.chemmap, {'Si': 1})
        self.assertEqual(value.nspecies, 1)

    def test_invalid_grain_metadata_refuses_before_deletion(self):
        for change in ['tag', 'centers', 'positions', 'cell']:
            value = sample([[4.9, 5, 5], [5.1, 5, 5]], [0, 1], pbc=True)
            if change == 'tag':
                value.polycrystal.set_tags([0, 99])
            elif change == 'centers':
                value.crystallitecenters[0, 0] = np.nan
            elif change == 'positions':
                value.polycrystal.positions[0, 0] = np.nan
            else:
                value.polycrystal.cell[0, 1] = 1
            with self.assertRaises(ValueError):
                prune(value, .5)
            self.assertEqual(len(value.polycrystal), 2)

    def test_generated_two_grain_configuration_prunes_and_writes(self):
        from grains import Grains
        from writefile import WriteData
        import tempfile
        numpy_state, random_state = np.random.get_state(), random.getstate()
        try:
            np.random.seed(17)
            random.seed(17)
            with contextlib.redirect_stdout(io.StringIO()):
                grains = Grains(np.diag([8., 8., 8.]), 2)
                value = grains.makepolycrystal(['Tungsten', 'Tungsten'], prunecriteria=0)
            initial = len(value.polycrystal)
            pairs = brute_pairs(value, .8)
            self.assertGreater(len(pairs), 0)
            prune(value, .8)
            self.assertLess(value.natoms, initial)
            self.assertEqual(brute_pairs(value, .8), [])
            with tempfile.TemporaryDirectory() as directory:
                filename = str(Path(directory) / 'tungsten')
                WriteData(filename).writelammpsdata(value)
                text = Path(filename + '.data').read_text()
                self.assertIn('%i atoms' % value.natoms, text)
                self.assertIn('1 atom types', text)
                rows = text.split('Atoms')[1].strip().splitlines()
                self.assertEqual(len(rows), value.natoms)
        finally:
            np.random.set_state(numpy_state)
            random.setstate(random_state)

    def test_seeded_geometries_match_brute_force_cross_grain_pairs(self):
        for pbc in [False, True, [True, False, True]]:
            for seed in range(8):
                rng = np.random.default_rng(seed)
                centers = rng.uniform(0, 10, (4, 3))
                positions = rng.uniform(0, 10, (70, 3))
                tags = np.argmin(np.sum((positions[:, None, :] - centers) ** 2, axis=2), axis=1)
                value = sample(positions, tags, centers, pbc=pbc)
                criteria = 1.7
                pairs = brute_pairs(value, criteria)
                boundary = set(value.grain_boundary_atoms(criteria))
                self.assertTrue(all(i in boundary and j in boundary for i, j in pairs))
                removed = set()
                for i, j in pairs:
                    if i not in removed and j not in removed:
                        removed.add(j)
                keep = [i for i in range(len(positions)) if i not in removed]
                prune(value, criteria)
                np.testing.assert_allclose(value.polycrystal.positions, positions[keep])
                self.assertEqual(brute_pairs(value, criteria), [])
                once = value.polycrystal.positions.copy()
                prune(value, criteria)
                np.testing.assert_array_equal(value.polycrystal.positions, once)


if __name__ == '__main__':
    unittest.main()
