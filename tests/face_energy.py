import itertools
import json
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from face_energy_obstruction import FiniteGroup, certificate, derive_fiber_group, exhaustive_face_equations


class FaceEnergyTests(unittest.TestCase):
    def test_all_minimal_d4_rules_against_independent_full_equations(self):
        group = derive_fiber_group()
        rules = json.loads((ROOT/'data/d4-pair-involutions.json').read_text())['rules']
        dimensions = set()
        for rule in rules:
            with self.subTest(rule=rule['id']):
                check = exhaustive_face_equations(group, rule['table'])
                dimensions.add(len(check['charge_basis']))
                self.assertEqual(check['states_checked'], 4096)
        self.assertEqual(dimensions, {1, 3, 4})

    def test_proper_nonabelian_quotient_and_nonzero_identity_index(self):
        elements = list(itertools.permutations(range(3)))
        group = FiniteGroup(elements[1:]+elements[:1])
        self.assertNotEqual(group.identity, 0)
        identity = list(range(group.n**2))
        self.assertEqual(certificate(group, identity)['normalized_face_charge_dimension'], 2)
        table = []
        for a, b in itertools.product(range(group.n), repeat=2):
            x = group.conj[b][a]
            y = group.mul[group.inv[x]][group.mul[a][b]]
            table.append(x*group.n+y)
        result = certificate(group, table)
        self.assertEqual(len(result['normal_subgroup']), 3)
        self.assertEqual(result['normalized_face_charge_dimension'], 1)
        self.assertEqual(len(exhaustive_face_equations(group, table)['charge_basis']), 1)

    def test_cyclic_four_nontrivial_update_fixes_quotient(self):
        group = FiniteGroup([tuple((i+a) % 4 for i in range(4)) for a in range(4)])
        table = [((a+2*b) % 4)*4+(-b) % 4 for a, b in itertools.product(range(4), repeat=2)]
        result = certificate(group, table)
        self.assertEqual(result['normal_subgroup'], [0, 2])
        self.assertEqual(result['element_charge_basis'], [[0, 1, 0, 1]])
        self.assertEqual(len(exhaustive_face_equations(group, table)['charge_basis']), 1)

    def test_random_product_preserving_maps_need_not_be_involutions(self):
        group = FiniteGroup([tuple((i+a) % 3 for i in range(3)) for a in range(3)])
        randomizer = random.Random(59378)
        for _ in range(12):
            table = []
            for a, b in itertools.product(range(3), repeat=2):
                x = randomizer.randrange(3)
                table.append(x*3+(a+b-x) % 3)
            exhaustive_face_equations(group, table)

    def test_invalid_rules_rejected(self):
        group = derive_fiber_group()
        table = list(range(64))
        table[0] = 1
        with self.assertRaisesRegex(ValueError, 'ordered product'):
            certificate(group, table)
        with self.assertRaisesRegex(ValueError, 'invalid entries'):
            certificate(group, table[:-1])

    def test_checked_census_excludes_nonconstant_faces_for_every_nonidentity_strict_rule(self):
        # Full table/braid revalidation and regeneration also run as a separate CI step.
        result = json.loads((ROOT/'data/d4-face-energy-obstruction.json').read_text())
        census = json.loads((ROOT/'data/d4-involution-braid-search.json').read_text())
        lookup = {s['id']: s for s in census['solutions'] if s['strict_braid']}
        ids = [i for row in result['strict_rule_classes'] for i in row['solution_ids']]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(set(ids), set(lookup))
        self.assertEqual(len(ids), 5731)
        for row in result['strict_rule_classes']:
            for i in row['solution_ids']:
                trivial = lookup[i]['table'] == list(range(64))
                self.assertEqual(row['normalized_face_charge_dimension'], 4 if trivial else 0)


if __name__ == '__main__':
    unittest.main()
