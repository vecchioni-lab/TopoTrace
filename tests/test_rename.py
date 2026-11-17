"""Renaming must preserve both public APIs and executable module aliases."""
import importlib
import inspect
import subprocess
import sys
import unittest


class RenameTests(unittest.TestCase):
    def test_new_and_legacy_modules_report_new_brand(self):
        for module in ('topotrace', 'topotrace.cgd', 'topostich', 'topostich.cgd'):
            with self.subTest(module=module):
                result = subprocess.run([sys.executable, '-m', module, '--version'],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), 'TopoTrace 0.5.0')

    def test_public_api_identity_and_packaged_catalog(self):
        for name, symbol in (('converter', 'convert_file'), ('converter', 'ConversionError'),
                             ('cgd', 'export_cgd'), ('nucleotides', 'recognize_residue'),
                             ('graph_viewer', 'marker_xml'), ('mmcif_input', 'read_mmcif')):
            with self.subTest(module=name, symbol=symbol):
                old = importlib.import_module('topostich.' + name)
                new = importlib.import_module('topotrace.' + name)
                self.assertIs(getattr(old, symbol), getattr(new, symbol))
        from topotrace.nucleotides import classify_atom
        self.assertEqual(classify_atom('DZ', 'C1').role, 'base')


class StitchingRenameTests(unittest.TestCase):
    def test_new_features_share_canonical_and_legacy_modules(self):
        for name, symbol in (('stitching', 'find_backbone_joins'), ('base_pairs', 'find_base_pairs'), ('asu_graphs', 'build_asu_graphs'), ('asu_view', 'render_html'), ('asu_png', 'render_png')):
            with self.subTest(module=name):
                canonical = importlib.import_module('topotrace.' + name)
                legacy = importlib.import_module('topostich.' + name)
                self.assertIs(getattr(legacy, symbol), getattr(canonical, symbol))

    def test_both_exporters_expose_stitching_and_graph_options(self):
        for package in ('topotrace', 'topostich'):
            with self.subTest(package=package):
                module = importlib.import_module(package + '.cgd')
                parameters = inspect.signature(module.export_cgd).parameters
                for name in ('stitch_nicks', 'stitch_terminal_stacks', 'asu_graphs',
                             'asu_chain_png', 'asu_residue_png', 'base_pairs'):
                    self.assertIn(name, parameters)
