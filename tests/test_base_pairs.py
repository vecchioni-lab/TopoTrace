"""C1'-prime base nodes and base-pair links (opt-in ``base_pairs``).

Each nucleotide with an O3' node gains a C1' base node joined to it. Base
pairs follow NAtopo's residue representation: hydrogen contacts (workbook
contact atoms, < 4.0 A, contact vector within 30 deg of both base planes,
different chain copies; plus deposited mmCIF ``hydrog`` links) grouped per
residue-pair symmetry orbit, and a pair is significant when each base
contributes at least two distinct atoms (``structural_features.FeatureGraph``).
Each significant pair becomes a C1'-C1' edge in its symmetry image.
"""
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from test_converter import crystal
from test_stitching import SCHEMAS, assert_reciprocal, nucleotide
from topotrace.cgd import export_cgd

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
RING = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')


def base(cx, cy, z, facing=1, *, dz=0.):
    """Cytosine in the plane z, N3 pointing along +x (facing=1) or -x (facing=-1).

    Exocyclic O2 (on C2) and N4 (on C4) and C1' (on N1) lie radially outside the
    ring. The backbone (P, O5', C5', C3', O3') sits behind C1'; dz offsets it
    along z so consecutive residues can be bonded O3'(i)-P(i+1).
    """
    def at(angle, radius, height=z):
        a = math.radians(angle)
        return (cx + facing*radius*math.cos(a), cy + radius*math.sin(a), height)
    atoms = {name: at(60*k - 120, 1.39) for k, name in enumerate(RING)}
    atoms.update(O2=at(-60, 2.62), N4=at(60, 2.62))
    atoms["C1'"] = at(-120, 2.86)
    bx, by = cx - facing*4.0, cy - 3.0
    for name, offset in (('P', -.8), ("O5'", -.2), ("C5'", .4), ("C3'", .2), ("O3'", .8)):
        atoms[name] = (bx, by, z + offset*dz if dz else z + offset)
    return atoms


def strand(serial, chain, cx, zs, facing, *, cy=10., drop=()):
    """Cytosines stacked at heights zs (5'->3'), consecutive O3'-P bonded at 1.8 A."""
    text = ''
    for i, z in enumerate(zs):
        d = 1 if zs[-1] >= zs[0] else -1
        atoms = base(cx, cy, z, facing, dz=d)
        if i == 0:
            del atoms['P']
        for name in drop:
            atoms.pop(name, None)
        text += nucleotide(serial + 20*i, i+1, chain, atoms, 'DC')
    return text + 'TER\n'


PAIR_GAP = 5.68  # ring centers; N3...N3 = 2.90 A, O2...O2 = N4...N4 = 3.06 A


def duplex(cell=40.):
    """Two antiparallel two-residue strands; A1 pairs B2 and A2 pairs B1."""
    return (crystal(cell, cell, cell) + strand(1, 'A', 10., (10., 13.4), 1)
            + strand(100, 'B', 10. + PAIR_GAP, (13.4, 10.), -1))


class BasePairTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.pdb'
        self.out = self.root / 'graph.cgd'

    def export(self, text, **kw):
        if text is not None:
            self.source.write_text(text)
        for path in [self.out, *(Path(p) for k, p in kw.items() if k not in ('base_pairs',) and isinstance(p, Path))]:
            if path.exists():
                path.unlink()
        kw.setdefault('base_pairs', True)
        return export_cgd(self.source, self.out, **kw)

    def pairs(self, r):
        return sorted((e['source'], e['target'], e['image_operation'])
                      for e in r['stitched_edges'] if e['kind'] == 'base_pair')

    def test_duplex_gains_c1_nodes_links_and_pairs(self):
        r = self.export(duplex())
        self.assertEqual(r['schema_version'], 3)
        self.assertEqual(r['base_atom'], "C1'")
        labels = [n['label'] for n in r['nodes']]
        self.assertEqual(labels, ['S1_N1', 'S1_N2', 'S2_N1', 'S2_N2', 'S1_B1', 'S1_B2', 'S2_B1', 'S2_B2'])
        self.assertEqual([n['node_type'] for n in r['nodes']], ['backbone']*4 + ['base']*4)
        self.assertEqual([n['atom'] for n in r['nodes'][4:]], ["C1'"]*4)
        self.assertEqual(r['base_links'], [['S1_N1', 'S1_B1'], ['S1_N2', 'S1_B2'], ['S2_N1', 'S2_B1'], ['S2_N2', 'S2_B2']])
        self.assertEqual(self.pairs(r), [('S1_B1', 'S2_B2', 'x,y,z'), ('S1_B2', 'S2_B1', 'x,y,z')])
        self.assertEqual([n['degree'] for n in r['nodes']], [2]*8)
        self.assertEqual(r['edge_count'], 2 + 4 + 2)
        self.assertTrue(all(n['base_paired'] for n in r['nodes'][4:]))
        summary = r['stitch_parameters']['base_pairs']
        self.assertEqual(summary['unpaired_base_nodes'], [])
        self.assertEqual(summary['hydrogen_cutoff_A'], 4.0)
        self.assertEqual(summary['min_atoms_per_base'], 2)
        statuses = sorted(c['status'] for c in r['base_pair_candidates'])
        self.assertEqual(statuses, ['stitched', 'stitched'])
        candidate = r['base_pair_candidates'][0]
        self.assertEqual(candidate['atoms_1'], ['N3', 'N4', 'O2'])
        self.assertEqual(candidate['atoms_2'], ['N3', 'N4', 'O2'])
        self.assertTrue(candidate['significant'])
        self.assertEqual(candidate['origin'], 'inferred')
        self.assertAlmostEqual(candidate['distance_min_A'], 2.9, places=2)
        edge = next(e for e in r['stitched_edges'] if e['kind'] == 'base_pair')
        self.assertAlmostEqual(edge['distance_A'], math.dist(r['nodes'][4]['cartesian'], r['nodes'][7]['cartesian']))
        assert_reciprocal(self, self.out.read_text(), 'P1')

    def test_stacked_neighbours_are_not_paired(self):
        # Cross-strand contacts between levels are > 4 A or > 30 deg off the planes.
        r = self.export(duplex())
        self.assertEqual(len(r['base_pair_candidates']), 2)

    def test_unpaired_base_is_kept_with_degree_one_and_flagged(self):
        lone = strand(200, 'C', 30., (30.,), 1)
        r = self.export(duplex() + lone)
        node = next(n for n in r['nodes'] if n['label'] == 'S3_B1')
        self.assertEqual((node['degree'], node['base_paired']), (1, False))
        unpaired = r['stitch_parameters']['base_pairs']['unpaired_base_nodes']
        self.assertEqual([u['label'] for u in unpaired], ['S3_B1'])
        self.assertTrue(any('unpaired' in w for w in r['warnings']))

    def test_one_atom_per_base_is_a_weak_contact(self):
        # Ring centers 6.7 A apart: only N3...N3 (3.92 A) is under 4 A.
        r = self.export(crystal(40, 40, 40) + strand(1, 'A', 10., (10.,), 1)
                        + strand(100, 'B', 16.7, (10.,), -1))
        candidate, = r['base_pair_candidates']
        self.assertEqual((candidate['status'], candidate['significant'], candidate['atom_pairs']),
                         ('weak_contact', False, 1))
        self.assertEqual(self.pairs(r), [])
        self.assertEqual(len(r['stitch_parameters']['base_pairs']['unpaired_base_nodes']), 2)
        r = self.export(None, min_atoms_per_base=1)
        self.assertEqual(self.pairs(r), [('S1_B1', 'S2_B1', 'x,y,z')])

    def test_contact_out_of_the_base_planes_is_rejected(self):
        # Raising B by 2.5 A tilts every contact ~40 deg out of both planes.
        text = (crystal(40, 40, 40) + strand(1, 'A', 10., (10.,), 1)
                + strand(100, 'B', 10. + PAIR_GAP, (12.5,), -1))
        r = self.export(text)
        self.assertEqual(r['base_pair_candidates'], [])
        r = self.export(None, hbond_plane_angle=45.)
        self.assertEqual(self.pairs(r), [('S1_B1', 'S2_B1', 'x,y,z')])
        r = self.export(None, hbond_plane_angle=None)
        self.assertEqual(len(self.pairs(r)), 1)

    def test_same_chain_copy_is_not_paired(self):
        # NAtopo's interchain rule: one strand folding back on itself in its own image.
        a = base(10., 10., 10., 1, dz=1)
        del a['P']
        b = base(10. + PAIR_GAP, 10., 10., -1, dz=-1)
        b['P'] = (b["O5'"][0], b["O5'"][1], 20.)
        text = crystal(40, 40, 40) + nucleotide(1, 1, 'A', a, 'DC') + nucleotide(30, 2, 'A', b, 'DC') + 'TER\n'
        r = self.export(text)
        self.assertEqual(r['base_pair_candidates'], [])

    def test_pair_with_own_twofold_mate_across_a_cell(self):
        # P2 (b axis): the base at x=12.16 pairs its mate -x+1,y,-z+1 at x=17.84.
        text = crystal(30, 30, 30, group='P 1 2 1') + strand(1, 'A', 12.16, (15.,), 1)
        r = self.export(text)
        edge, = [e for e in r['stitched_edges'] if e['kind'] == 'base_pair']
        self.assertEqual((edge['source'], edge['target'], edge['image_operation']), ('S1_B1', 'S1_B1', '-x+1,y,-z+1'))
        node = next(n for n in r['nodes'] if n['label'] == 'S1_B1')
        self.assertEqual(node['degree'], 2)  # its O3' link and one self-image neighbour
        candidate, = r['base_pair_candidates']
        self.assertEqual(candidate['atoms_1'], candidate['atoms_2'])
        assert_reciprocal(self, self.out.read_text(), 'P2')

    def test_residue_without_c1_has_no_base_node(self):
        text = (crystal(40, 40, 40) + strand(1, 'A', 10., (10.,), 1, drop=("C1'",))
                + strand(100, 'B', 10. + PAIR_GAP, (10.,), -1))
        r = self.export(text)
        self.assertEqual([n['label'] for n in r['nodes']], ['S1_N1', 'S2_N1', 'S2_B1'])
        candidate, = r['base_pair_candidates']
        self.assertEqual(candidate['status'], 'no_graph_node')
        missing = r['stitch_parameters']['base_pairs']['residues_without_base_node']
        self.assertEqual([m['seq'] for m in missing], ['1'])

    def test_deposited_hydrogen_link_counts(self):
        # Geometric screen finds only N3...N3; a deposited O2...N4 hydrog link adds one atom per base.
        rows, serial = [], 0
        for label, cx, facing in (('A', 10., 1), ('B', 16.7, -1)):
            atoms = base(cx, 10., 10., facing, dz=1)
            del atoms['P']
            for name, xyz in atoms.items():
                serial += 1
                rows.append(f'ATOM {serial} {name[0]} "{name}" . DC {label} 1 {label} 1 ? '
                            f'{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} 1 20 1\n')
        from test_mmcif_input import HEADER, TAGS
        conn = ('loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n'
                '_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n'
                '_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n'
                'h1 hydrog A 1 O2 B 1 N4\n')
        self.source = self.root / 'input.cif'
        header = HEADER.replace('_cell.length_a 20', '_cell.length_a 40').replace(
            '_cell.length_b 30', '_cell.length_b 40')
        self.source.write_text(header + TAGS + ''.join(rows) + conn)
        r = self.export(None)
        candidate, = r['base_pair_candidates']
        self.assertEqual((candidate['status'], candidate['origin'], candidate['atom_pairs']), ('stitched', 'both', 2))
        self.assertEqual((candidate['atoms_1'], candidate['atoms_2']), (['N3', 'O2'], ['N3', 'N4']))

    def test_base_pairs_with_nick_and_stack_stitching(self):
        r = self.export(duplex(), stitch_nicks=True, stitch_terminal_stacks=True)
        self.assertEqual(r['schema_version'], 3)
        self.assertEqual(r['stitched_edge_counts'], {'base_pair': 2})
        # Terminal stacks stay on O3' nodes: no stack between these ends exists here.
        self.assertTrue(all(c['source'] is None or '_N' in c['source'] for c in r['terminal_stack_candidates']))

    def test_default_and_stitched_outputs_have_no_base_nodes(self):
        r = self.export(duplex(), base_pairs=False)
        self.assertEqual(r['schema_version'], 1)
        self.assertNotIn('base_links', r)
        self.assertEqual(len(r['nodes']), 4)
        r = self.export(None, base_pairs=False, stitch_nicks=True)
        self.assertEqual(r['schema_version'], 2)
        self.assertNotIn('node_type', r['nodes'][0])

    def test_views_and_asu_graphs(self):
        cmm, crystal_cmm, asu = self.root/'a.cmm', self.root/'s.cmm', self.root/'asu.json'
        r = self.export(duplex(), cmm=cmm, crystal_cmm=crystal_cmm, asu_graphs=asu)
        root = ET.fromstring(cmm.read_text())
        names = [m.get('name') for m in root.iter('marker_set')]
        self.assertIn('input_O3_C1_base_paired base pairs', names)
        segment = next(m for m in root.iter('marker_set') if m.get('name').endswith('segment 1'))
        self.assertEqual(len(segment.findall('marker')), 4)  # two O3' and two C1' markers
        self.assertEqual(len(segment.findall('link')), 3)    # one backbone link and two O3'-C1' links
        doc = json.loads(asu.read_text())
        residue = doc['residue_graph']
        self.assertEqual(residue['edge_type_counts'], {'base_pair': 2, 'covalent': 2})
        self.assertEqual(doc['chain_graph']['edge_type_counts'], {'base_pair': 1})
        self.assertEqual(doc['stitching']['base_pairs'], True)
        try:
            import jsonschema
        except ImportError:
            self.skipTest('jsonschema not installed')
        jsonschema.validate(doc, json.loads((SCHEMAS/'asu-graphs.schema.json').read_text()))

    def test_cli_and_receipt_schema_3(self):
        self.source.write_text(duplex())
        receipt = self.root/'graph.json'
        run = subprocess.run([sys.executable, '-m', 'topotrace.cgd', str(self.source), '--out', str(self.out),
                              '--base-pairs', '--receipt', str(receipt)], text=True, capture_output=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn('2 base pairs', run.stdout)
        self.assertIn("4 C1'", run.stdout)
        r = json.loads(receipt.read_text())
        self.assertEqual(r['schema_version'], 3)
        try:
            import jsonschema
        except ImportError:
            self.skipTest('jsonschema not installed')
        jsonschema.validate(r, json.loads((SCHEMAS/'cgd-receipt-3.schema.json').read_text()))

    def test_systre_audit_counts_base_links(self):
        spec = importlib.util.spec_from_file_location('check_stitched_cgd', TOOLS/'check_stitched_cgd.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        r = self.export(duplex())
        full, half = module.expected_pairs(r)
        self.assertEqual(sum(full.values()), 2 + 4 + 2)
        self.assertEqual(full[frozenset(('S1_N1', 'S1_B1'))], 1)


if __name__ == '__main__':
    unittest.main()
