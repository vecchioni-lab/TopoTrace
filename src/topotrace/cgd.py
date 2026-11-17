"""Chain-aware O3-prime graphs in the named-source CRYSTAL CGD dialect.

This is a separate exporter: it does not change the existing CIF/PDB/mmCIF
converter. Explicit edges represent traces, not chemical bond assertions.
Optional nick (``stitch_nicks``) and chain-end stack (``stitch_terminal_stacks``)
stitching add symmetry-aware joins found by the rules in ``stitching.py``; the
default is unchanged. Optional ``base_pairs`` adds a C1-prime base node per
nucleotide, linked to its O3-prime node, and joins C1-prime nodes of NAtopo
significant hydrogen pairs (``base_pairs.py``). Optional ASU quotient graphs (chain and residue
resolution, covalent bonds and stitched joins) come from ``asu_graphs.py``.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import gemmi

from . import __version__
from .converter import ConversionError, _fractional_transform, _select, read_structure
from .asu_graphs import build_asu_graphs, build_atom_graph, numbered_operations
from .asu_view import render_html
from .graph_viewer import marker_xml
from .stitching import IDENTITY_OP, StitchContext, find_backbone_joins, find_terminal_stacks
from .base_pairs import BASE_ATOM, base_nodes, find_base_pairs


def _space_group(pdb):
    group = gemmi.find_spacegroup_by_name(pdb.space_group)
    if group is None:
        raise ConversionError(f"Unknown space group: {pdb.space_group}")
    if not gemmi.UnitCell(*pdb.cell).is_compatible_with_spacegroup(group):
        raise ConversionError("Cell is incompatible with the named space-group setting")
    operations = list(group.operations())
    if pdb.symops:
        try:
            supplied = {gemmi.Op(op).wrap().triplet() for op in pdb.symops.values()}
        except (RuntimeError, ValueError) as error:
            raise ConversionError(f"Invalid input symmetry operation: {error}") from error
        expected = {op.wrap().triplet() for op in operations}
        if supplied != expected:
            raise ConversionError("Input symmetry operations disagree with the complete named space-group setting")
    return group, operations


def _check_positions(fracs, operations):
    # Systre compares points modulo its primitive lattice at 0.001 fractional
    # precision. This conservative conventional-cell guard (0.003 in each
    # coordinate) also catches near-coincidences after centering transformations.
    # Special positions need stabilizer-aware edge degrees: defer them explicitly.
    tolerance = .003
    identity = gemmi.Op('x,y,z').triplet()
    for i, point in enumerate(fracs):
        for op in operations:
            image = op.apply_to_xyz(point)
            for j in range(i + 1):
                if i == j and op.triplet() == identity:
                    continue
                delta = [a-b for a,b in zip(image,fracs[j])]
                if max(abs(d-round(d)) for d in delta) <= tolerance:
                    kind = 'special position' if i == j else 'symmetry-equivalent or near-coincident nodes'
                    raise ConversionError(f"CGD cannot preserve separate trace identities at {kind} "
                                          f"(nodes {j+1}, {i+1}; operation {op.triplet()})")


def _number(value):
    return format(value, '.12g')


def _neighbor_rows(nodes, edges, stitched=()):
    """CGD neighbor-image rows per node label, without repeated rows.

    Periodic joins use the inverse affine image for the reciprocal row, not a
    copied or negated translation. A join from a node to its own image under an
    involution (e.g. a base stacked on its 2-fold mate) has identical forward
    and reciprocal rows: it is one neighbor, written once and counted once.
    """
    neighbors = {node['label']: [] for node in nodes}
    seen = {node['label']: set() for node in nodes}
    by_label = {node['label']: node for node in nodes}

    def add(label, image):
        key = tuple(round(x, 9) for x in image)
        if key not in seen[label]:
            seen[label].add(key)
            neighbors[label].append(image)
    for left, right in edges:
        add(left, by_label[right]['fractional'])
        add(right, by_label[left]['fractional'])
    for edge in stitched:
        add(edge['source'], edge['target_image_fractional'])
        add(edge['target'], edge['reverse_image_fractional'])
    return neighbors


def _serialize(graph):
    lines = ['crystal', 'name '+graph['name'],
             'cell '+' '.join(_number(x) for x in graph['cell']),
             'group '+graph['space_group']]
    neighbors = _neighbor_rows(graph['nodes'], graph['edges'] + graph.get('base_links', []),
                               graph.get('stitched_edges', ()))
    for node in graph['nodes']:
        label = node['label']
        lines.append(f"atom {label} {node['degree']} "+' '.join(map(_number,node['fractional'])))
        for image in neighbors[label]:
            lines.append(f"edge {label} "+' '.join(map(_number,image)))
    lines.append('end')
    return '\n'.join(lines)+'\n'


def _triplet(g):
    from .asu_graphs import to_gemmi
    return to_gemmi(g).triplet()


def _write_new_files(payloads):
    paths = [path.resolve() for path in payloads]
    if len(paths) != len(set(paths)):
        raise ConversionError('Output paths must be distinct')
    for path in payloads:
        if path.exists() or path.is_symlink():
            raise ConversionError(f'Refusing to overwrite output: {path}')
        if not path.parent.is_dir():
            raise ConversionError(f'Output parent does not exist: {path.parent}')
    created = []
    try:
        for path, content in payloads.items():
            with path.open('xb') as stream:
                created.append(path)
                stream.write(content)
    except OSError as error:
        for path in created:
            path.unlink()
        raise ConversionError(str(error)) from error


def export_cgd(source, output, *, cmm=None, crystal_cmm=None, receipt=None,
               name=None, model=1, altloc=None, input_format='auto',
               stitch_nicks=False, join_cutoff=3.5, o5_cutoff=6.0, include_o5_fallback=False,
               stitch_terminal_stacks=False, asu_graphs=None, asu_chain_html=None, asu_residue_html=None,
               asu_chain_png=None, asu_residue_png=None, base_pairs=False, hbond_cutoff=4.0,
               hbond_plane_angle=30.0, min_atoms_per_base=2, asu_atom_html=None, asu_atom_png=None):
    """Write a new CGD and optional ChimeraX views/JSON receipt; return receipt.

    Every output is self-contained and refused if it exists. cmm contains the
    observed asymmetric-unit trace. crystal_cmm contains one intact copy for
    each conventional-cell symmetry operation.

    With stitch_nicks, unambiguous P...O3' backbone nick joins between segment
    termini (NAtopo rule, any symmetry image) become periodic edges. With
    stitch_terminal_stacks, pi stacks between chain ends (any polarity, any
    image) do too, unless the two ends are already backbone joined. Either
    option makes the receipt schema 2. Other candidates are reported, never
    silently used.

    With base_pairs, every nucleotide with an O3' node also gets a C1' base
    node joined to it, and significant hydrogen pairs (NAtopo residue
    representation: hydrogen contacts < hbond_cutoff within hbond_plane_angle of
    both base planes, at least min_atoms_per_base distinct atoms per base)
    become C1'-C1' edges in any symmetry image. The receipt becomes schema 3.
    Unpaired bases keep a degree-1 C1' node and are flagged.

    asu_graphs (JSON), asu_chain_html and asu_residue_html write the ASU
    quotient graphs of the covalent bonds and stitched joins; asu_chain_png and
    asu_residue_png write PNG images of them (needs matplotlib). asu_atom_html
    and asu_atom_png draw a third view at atom resolution: the CGD's own O3'
    and C1' nodes with every trace edge, O3'-C1' link, stitched join and base pair.
    """
    stitching = bool(stitch_nicks or stitch_terminal_stacks or base_pairs)
    asu_requested = any(p is not None for p in (asu_graphs, asu_chain_html, asu_residue_html,
                                                 asu_chain_png, asu_residue_png, asu_atom_html, asu_atom_png))
    source, output = Path(source), Path(output)
    targets = [Path(p) for p in (output,cmm,crystal_cmm,receipt,asu_graphs,asu_chain_html,asu_residue_html,
                                 asu_chain_png,asu_residue_png,asu_atom_html,asu_atom_png)
               if p is not None]
    # Check before dict construction so identical companion paths cannot collapse.
    if len({p.resolve() for p in targets}) != len(targets):
        raise ConversionError('Output paths must be distinct')
    if any(p.resolve() == source.resolve() for p in targets):
        raise ConversionError('Input cannot be an output')
    if name is None:
        suffix = {(False, False): '_O3_chain_aware', (True, False): '_O3_nick_stitched',
                  (False, True): '_O3_end_stack_stitched', (True, True): '_O3_stitched'}
        tail = suffix[bool(stitch_nicks), bool(stitch_terminal_stacks)]
        if base_pairs:
            tail = '_O3_C1_base_paired' if tail == '_O3_chain_aware' else tail.replace('_O3', '_O3_C1', 1)
        name = re.sub(r'[^A-Za-z0-9_.-]', '_', source.stem)+tail
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*',name):
        raise ConversionError('Graph name must be a single ASCII identifier')
    raw = source.read_bytes()
    pdb = read_structure(raw.decode('utf-8'),model=model,input_format=input_format)
    atoms, pairs, _, missing = _select(pdb,["O3'"],altloc)
    transform = _fractional_transform(pdb)
    group, operations = _space_group(pdb)
    fracs = [[sum(row[k]*atom.xyz[k] for k in range(3))+row[3] for row in transform]
             for atom in atoms]
    _check_positions(fracs,operations)
    counters, labels = Counter(), []
    for atom in atoms:
        counters[atom.segment] += 1
        labels.append(f'S{atom.segment}_N{counters[atom.segment]}')
    stitched, candidates, stacks, nick_summary, stack_summary = [], [], [], None, None
    context = None
    if stitching or asu_requested:
        context = StitchContext(pdb, altloc=altloc, node_atoms=atoms, node_labels=labels,
                                transform=transform, operations=operations, trace_pairs=pairs)
        if stitch_nicks:
            edges, candidates, nick_summary = find_backbone_joins(
                context, join_cutoff=join_cutoff, o5_cutoff=o5_cutoff,
                include_o5_fallback=include_o5_fallback)
            stitched += edges
        if stitch_terminal_stacks:
            edges, stacks, stack_summary = find_terminal_stacks(context)
            stitched += edges
    node_atoms, node_labels, node_fracs = list(atoms), list(labels), list(fracs)
    links, base_missing, pair_records, pair_summary = [], [], [], None
    if base_pairs:
        # C1' base nodes follow all O3' nodes; S<segment>_B<k> belongs to S<segment>_N<k>.
        found, base_missing = base_nodes(context, atoms)
        base_node_of = {}
        for i, c1 in found:
            index = len(context.node_atoms)
            context.node_atoms.append(c1)
            context.node_labels.append(labels[i].replace('_N', '_B', 1))
            base_node_of[c1.residue_key] = index
            links.append((i, index))
            context._occupy(i, index, IDENTITY_OP)
        node_atoms, node_labels = context.node_atoms, context.node_labels
        node_fracs = fracs + [[sum(row[k]*c1.xyz[k] for k in range(3))+row[3] for row in transform]
                              for _, c1 in found]
        _check_positions(node_fracs, operations)
        edges, pair_records, pair_summary = find_base_pairs(
            context, base_node_of, hydrogen_cutoff=hbond_cutoff, max_plane_angle=hbond_plane_angle,
            min_atoms_per_base=min_atoms_per_base)
        stitched += edges
    base_links = [[node_labels[a], node_labels[b]] for a, b in links]
    rows = _neighbor_rows([dict(label=label, fractional=frac) for label, frac in zip(node_labels, node_fracs)],
                          [[labels[a], labels[b]] for a, b in pairs] + base_links, stitched)
    nodes = [dict(label=node_labels[i],degree=len(rows[node_labels[i]]),segment=a.segment,
                  chain=a.chain,seq=a.seq,insertion=a.insertion,residue=a.residue,
                  atom=a.name,altloc=a.altloc,label_asym=a.label_asym,
                  label_seq=a.label_seq,label_comp=a.label_comp,label_atom=a.label_atom,
                  cartesian=list(a.xyz),fractional=node_fracs[i],description=a.description)
             for i,a in enumerate(node_atoms)]
    unpaired = []
    if base_pairs:
        paired = {label for e in stitched if e['kind'] == 'base_pair' for label in (e['source'], e['target'])}
        base_labels = {b for _, b in base_links}
        for node in nodes:
            node['node_type'] = 'base' if node['label'] in base_labels else 'backbone'
            if node['node_type'] == 'base':
                node['base_paired'] = node['label'] in paired
                if not node['base_paired']:
                    unpaired.append(dict(label=node['label'], description=node['description']))
        pair_summary.update(base_node_count=len(links), residues_without_base_node=base_missing,
                            unpaired_base_nodes=unpaired)
    if stitching:
        kinds = Counter(e['kind'] for e in stitched)
        warnings = []
        if base_pairs:
            counts = Counter(c['status'] for c in pair_records)
            warnings.append(f"Base pairs: {kinds['base_pair']} C1'-C1' links from significant hydrogen pairs; "
                            f"{len(links)} C1' base nodes linked to their O3' nodes.")
            for status, text in (('weak_contact', 'hydrogen contacts with fewer than min_atoms_per_base atoms per base were not linked'),
                                 ('no_graph_node', "base pairs without a C1' or O3' graph node were not linked"),
                                 ('duplicate_edge', 'base pairs duplicating an existing edge were not linked')):
                if counts[status]:
                    warnings.append(f'{counts[status]} {text}.')
            if unpaired:
                warnings.append(f"{len(unpaired)} unpaired C1' base nodes kept (flagged base_paired=false; "
                                "degree 1 unless otherwise linked).")
            if base_missing:
                warnings.append(f"{len(base_missing)} residues with an O3' node lack C1' and have no base node.")
        if stitch_nicks:
            counts = Counter(c['status'] for c in candidates)
            nick_count = len(stitched) - kinds['terminal_stack'] - kinds['base_pair']
            warnings.append(f'Nick stitching: {nick_count} backbone joins stitched.')
            for status, text in (('ambiguous', 'ambiguous nick candidates (several partners) were not stitched'),
                                 ('segment_internal', 'nick candidates at segment-internal breaks were not stitched'),
                                 ('o5_fallback_disabled', "O5' fallback nick candidates reported only (use include_o5_fallback)"),
                                 ('no_graph_node', 'nick candidates without an O3-prime graph node were not stitched'),
                                 ('duplicate_edge', 'nick candidates duplicating an existing edge were not stitched')):
                if counts[status]:
                    warnings.append(f'{counts[status]} {text}.')
        if stitch_terminal_stacks:
            counts = Counter(c['status'] for c in stacks)
            same = sum(e['evidence']['same_polarity'] for e in stitched if e['kind'] == 'terminal_stack')
            warnings.append(f"Chain-end stack stitching: {kinds['terminal_stack']} terminal stacks stitched"
                            + (f' ({same} same-polarity 5\'-5\' or 3\'-3\').' if same else '.'))
            for status, text in (('backbone_joined', 'chain-end stacks between ends already backbone joined were not stitched'),
                                 ('no_graph_node', 'chain-end stacks without an O3-prime graph node were not stitched'),
                                 ('duplicate_edge', 'chain-end stacks duplicating an existing edge were not stitched')):
                if counts[status]:
                    warnings.append(f'{counts[status]} {text}.')
            issues = stack_summary['plane_diagnostics']
            if issues:
                warnings.append(f'{len(issues)} chain-end residues lack a usable base-ring plane; see stitch_parameters.')
        if not stitch_terminal_stacks:
            warnings.append('Chain-end pi stacks are not stitched.')
        if not stitch_nicks:
            warnings.append('Backbone nicks are not stitched.')
    else:
        warnings = ['Unstitched finite strand components; no nick, junction, or symmetry-image joins.']
    if missing:
        warnings.append(f'{len(missing)} missing O3-prime selections skipped; retained nodes bridge within each segment.')
    if not pairs and not stitched and not links:
        warnings.append('No edges: all declared degrees are zero; Systre may report its nearest-node fallback.')
    graph = dict(schema_version=3 if base_pairs else 2 if stitching else 1,format='Systre CGD CRYSTAL',tool_version=__version__,
                 name=name,source_name=source.name,source_sha256=hashlib.sha256(raw).hexdigest(),
                 source_format=pdb.source_format,model=model,altloc=altloc,representative_atom="O3'",
                 cell=list(pdb.cell),space_group=group.xhm().replace(' ',''),
                 symmetry_operations=[op.triplet() for op in operations],
                 node_count=len(nodes),edge_count=len(pairs)+len(links)+len(stitched),
                 directed_edge_rows=2*(len(pairs)+len(links)+len(stitched)),
                 segment_count=len(counters),nodes=nodes,
                 edges=[[labels[a],labels[b]] for a,b in pairs],missing_atoms=missing,
                 stitching=stitching,warnings=warnings)
    if stitching:
        graph.update(trace_edge_count=len(pairs),stitched_edge_count=len(stitched),
                     stitched_edge_counts=dict(sorted(Counter(e['kind'] for e in stitched).items())),
                     stitched_edges=stitched,join_candidates=candidates,terminal_stack_candidates=stacks,
                     stitch_parameters=dict(nicks=nick_summary,terminal_stacks=stack_summary))
    if base_pairs:
        graph.update(base_atom=BASE_ATOM,base_node_count=len(links),base_link_count=len(links),
                     base_links=base_links,base_pair_candidates=pair_records)
        graph['stitch_parameters']['base_pairs'] = pair_summary
    payloads = {output:_serialize(graph).encode()}
    if asu_requested:
        quotients = build_asu_graphs(context, numbered_operations(pdb, operations), name=name,
                                     space_group=graph['space_group'])
        document = dict(schema_version=1, format='TopoTrace ASU quotient graphs', tool_version=__version__,
                        name=name, source_name=source.name, source_sha256=graph['source_sha256'],
                        model=model, altloc=altloc, cell=list(pdb.cell), space_group=graph['space_group'],
                        symmetry_operations=[dict(id=i, operation=_triplet(g))
                                             for i, g in numbered_operations(pdb, operations)],
                        stitching=dict(nicks=bool(stitch_nicks), terminal_stacks=bool(stitch_terminal_stacks),
                                       **(dict(base_pairs=True) if base_pairs else {})),
                        chain_graph=quotients['chain'], residue_graph=quotients['residue'])
        if asu_graphs is not None:
            payloads[Path(asu_graphs)] = (json.dumps(document, indent=2, allow_nan=False)+'\n').encode()
        if asu_chain_html is not None:
            payloads[Path(asu_chain_html)] = render_html(quotients['chain']).encode()
        if asu_residue_html is not None:
            payloads[Path(asu_residue_html)] = render_html(quotients['residue']).encode()
        if asu_chain_png is not None or asu_residue_png is not None:
            from .asu_png import render_png
            for path, resolution in ((asu_chain_png, 'chain'), (asu_residue_png, 'residue')):
                if path is not None:
                    payloads[Path(path)] = render_png(quotients[resolution])
        if asu_atom_html is not None or asu_atom_png is not None:
            quotients['atom'] = build_atom_graph(graph, quotients['residue'], numbered_operations(pdb, operations),
                                                 name=name, space_group=graph['space_group'])
            if asu_atom_html is not None:
                payloads[Path(asu_atom_html)] = render_html(quotients['atom']).encode()
            if asu_atom_png is not None:
                from .asu_png import render_png
                payloads[Path(asu_atom_png)] = render_png(quotients['atom'])
        graph['asu_graphs'] = {resolution: dict(nodes=len(q['nodes']), edges=len(q['edges']),
                                                edge_type_counts=q['edge_type_counts'],
                                                periodicity=q['summary'])
                               for resolution, q in quotients.items()}
    if cmm is not None:
        payloads[Path(cmm)] = marker_xml(graph).encode()
    if crystal_cmm is not None:
        payloads[Path(crystal_cmm)] = marker_xml(graph,crystal=True).encode()
    graph['artifacts'] = [dict(name=p.name,sha256=hashlib.sha256(data).hexdigest())
                          for p,data in payloads.items()]
    if receipt is not None:
        payloads[Path(receipt)] = (json.dumps(graph,indent=2,allow_nan=False)+'\n').encode()
    _write_new_files(payloads)
    return graph


def main(argv=None):
    parser = argparse.ArgumentParser(prog='topotrace-cgd', description='O3-prime trace to Systre CRYSTAL CGD; chains stay separate unless --stitch-nicks/--stitch-terminal-stacks/--base-pairs.')
    parser.add_argument('input',type=Path,help='PDB or native PDBx/mmCIF')
    parser.add_argument('--out','-out',required=True,type=Path,help='new .cgd file')
    parser.add_argument('--cmm',type=Path,help='ChimeraX asymmetric-unit graph view')
    parser.add_argument('--crystal-cmm',type=Path,help='ChimeraX intact symmetry copies (plus stitched joins)')
    parser.add_argument('--receipt',type=Path,help='save node identities, edges and provenance as JSON')
    parser.add_argument('--name',help='CGD identifier; default derived from source basename')
    parser.add_argument('--input-format',choices=('auto','pdb','mmcif'),default='auto')
    parser.add_argument('--model',type=int,default=1)
    parser.add_argument('--altloc')
    parser.add_argument('--stitch-nicks',action='store_true',
                        help="join unambiguous loose P...O3' backbone nicks across chains/symmetry images (NAtopo rule)")
    parser.add_argument('--join-cutoff',type=float,default=3.5,help="P...O3' cutoff in Angstrom (default 3.5)")
    parser.add_argument('--o5-cutoff',type=float,default=6.0,help="missing-phosphate O5'...O3' cutoff (default 6.0)")
    parser.add_argument('--include-o5-fallback',action='store_true',
                        help="also stitch unambiguous O5' fallback joins (reported only by default)")
    parser.add_argument('--stitch-terminal-stacks',action='store_true',
                        help="join chain ends whose bases pi-stack (5'-3', 5'-5' or 3'-3'; any symmetry image) unless already backbone joined")
    parser.add_argument('--base-pairs',action='store_true',
                        help="add a C1' base node per nucleotide (linked to its O3' node) and join C1' nodes of "
                             "significant hydrogen pairs (NAtopo rule; any symmetry image)")
    parser.add_argument('--hbond-cutoff',type=float,default=4.0,
                        help='base-pair hydrogen contact cutoff in Angstrom (default 4.0)')
    parser.add_argument('--hbond-plane-angle',type=float,default=30.0,
                        help='maximum contact-vector angle to each base plane in degrees (default 30; 90 disables)')
    parser.add_argument('--min-pair-atoms',type=int,default=2,
                        help='distinct contact atoms required on each base for a base pair (default 2)')
    parser.add_argument('--asu-graphs',type=Path,help='save the ASU quotient graphs (chain and residue resolution) as JSON')
    parser.add_argument('--asu-chain-html',type=Path,help='interactive HTML view of the chain-resolution ASU graph')
    parser.add_argument('--asu-residue-html',type=Path,help='interactive HTML view of the residue-resolution ASU graph')
    parser.add_argument('--asu-chain-png',type=Path,help='PNG image of the chain-resolution ASU graph (needs matplotlib)')
    parser.add_argument('--asu-residue-png',type=Path,help='PNG image of the residue-resolution ASU graph (needs matplotlib)')
    parser.add_argument('--asu-atom-html',type=Path,
                        help="interactive HTML view at atom resolution: O3' and C1' nodes and all their links")
    parser.add_argument('--asu-atom-png',type=Path,
                        help="PNG image of the atom-resolution view (needs matplotlib)")
    parser.add_argument('--json',action='store_true',help='print the receipt as JSON')
    parser.add_argument('--version',action='version',version=f'TopoTrace {__version__}')
    args=parser.parse_args(argv)
    try:
        result=export_cgd(args.input,args.out,cmm=args.cmm,crystal_cmm=args.crystal_cmm,
                          receipt=args.receipt,name=args.name,model=args.model,
                          altloc=args.altloc,input_format=args.input_format,
                          stitch_nicks=args.stitch_nicks,join_cutoff=args.join_cutoff,
                          o5_cutoff=args.o5_cutoff,include_o5_fallback=args.include_o5_fallback,
                          stitch_terminal_stacks=args.stitch_terminal_stacks,asu_graphs=args.asu_graphs,
                          asu_chain_html=args.asu_chain_html,asu_residue_html=args.asu_residue_html,
                          asu_chain_png=args.asu_chain_png,asu_residue_png=args.asu_residue_png,
                          base_pairs=args.base_pairs,hbond_cutoff=args.hbond_cutoff,
                          hbond_plane_angle=args.hbond_plane_angle,min_atoms_per_base=args.min_pair_atoms,
                          asu_atom_html=args.asu_atom_html,asu_atom_png=args.asu_atom_png)
    except (ConversionError,OSError,UnicodeError) as error:
        parser.exit(2,f'topotrace-cgd: {error}\n')
    if args.json:
        print(json.dumps(result,indent=2,allow_nan=False))
    elif result['stitching']:
        kinds = result['stitched_edge_counts']
        nicks = result['stitched_edge_count'] - kinds.get('terminal_stack', 0) - kinds.get('base_pair', 0)
        parts = ([f"{nicks} stitched nick joins"] if args.stitch_nicks else []) + \
                ([f"{kinds.get('terminal_stack', 0)} stitched chain-end stacks"] if args.stitch_terminal_stacks else [])
        nodes = f"{result['node_count']} nodes"
        if args.base_pairs:
            bases = result['base_node_count']
            nodes += f" ({result['node_count'] - bases} O3', {bases} C1')"
            parts = [f"{result['base_link_count']} O3'-C1'"] + parts + [f"{kinds.get('base_pair', 0)} base pairs"]
        print(f"Wrote {args.out}: {nodes}, {result['edge_count']} edges "
              f"({result['trace_edge_count']} trace, {', '.join(parts)}), "
              f"{result['segment_count']} strand segments")
    else:
        print(f"Wrote {args.out}: {result['node_count']} nodes, {result['edge_count']} edges, "
              f"{result['segment_count']} separate strand segments")
    if not args.json:
        for resolution, quotient in result.get('asu_graphs', {}).items():
            print(f"ASU {resolution} graph: {quotient['nodes']} nodes, {quotient['edges']} edge orbits")
            for line in quotient['periodicity']:
                print('  '+line)
        for warning in result['warnings']:
            print('Note: '+warning,file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
