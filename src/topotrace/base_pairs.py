"""C1'-prime base nodes and base-pair links for the O3-prime CGD graph.

Opt-in (``export_cgd(..., base_pairs=True)``, ``topotrace-cgd --base-pairs``).
Each nucleotide that has an O3' node gains a second node at its C1' atom,
joined to its own O3' node. Base pairs then join C1' nodes, possibly in another
symmetry image. The rule is NAtopo's residue representation, ported without
NetworkX, NumPy, SciPy or openpyxl:

* **Hydrogen contacts** (``nucleotide_interactions.find_hydrogen_bond_candidates``
  with the active ``nucleotide_interactions.json``): contact atoms from the lab's
  ``bases.xlsx`` Atom ID lists (A: N1 N3 N6 N7; G: N1 N2 N3 O6 N7; C: O2 N3 N4;
  T/U: O2 N3 O4; Donor/Acceptor columns ignored) on *different chain copies*
  (another segment, or any other symmetry/lattice image) closer than
  ``hydrogen_cutoff`` (4.0 A), with the contact vector within
  ``max_plane_angle`` (30 deg) of both base planes. Pairs already covalently
  linked are skipped. Deposited mmCIF ``_struct_conn`` ``hydrog`` links between
  two base-bearing residues count too (``crystal_interactions``: deposited and
  inferred layers count together, duplicates once).
* **Significant pair** (``structural_features.FeatureGraph._hydrogen``): contacts
  are grouped per residue-pair symmetry orbit; the pair is significant when each
  base contributes at least ``min_atoms_per_base`` (2) distinct atoms. A
  residue paired with its own two-fold copy counts each stored atom pair x...y'
  also as y...x'. Other orbits are reported as ``weak_contact``.

Every significant pair becomes a C1'-C1' edge in its symmetry image. Bases with
no significant pair keep their C1' node (degree 1) and are flagged unpaired.
Ring atom names, contact names of modified residues and base planes come from
TopoTrace's NARestraints-derived residue catalog (see ``stitching.base_planes``).
Every symmetry operation and lattice translation within reach is enumerated
exactly. No atoms are moved, added or repaired; links are topological, not
chemical bonds.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import math

from .converter import ConversionError, normalized_name
from .stitching import (IDENTITY, IDENTITY_OP, _image_op, _norm, _positive, _residue_identity, _shift,
                        _transform_plane, _wrapped, base_planes, op_key)

BASE_PAIR_METHOD = 'natopo_significant_hydrogen_pair_v1'
BASE_ATOM = "C1'"
# NAtopo ../manuscripts and spreadsheet/bases.xlsx, Atom ID column of each sheet.
CONTACT_ATOMS = {'A': ('N1', 'N3', 'N6', 'N7'), 'G': ('N1', 'N2', 'N3', 'O6', 'N7'),
                 'C': ('O2', 'N3', 'N4'), 'T': ('O2', 'N3', 'O4'), 'U': ('O2', 'N3', 'O4')}
CONTACT_ELEMENTS = {'N', 'O', 'S', 'SE'}
DEFAULTS = dict(hydrogen_cutoff_A=4.0, max_plane_angle_deg=30.0, min_atoms_per_base=2)


def _contact_names(code, observed):
    """(canonical -> actual contact atom names, base letter) or (None, issue)."""
    from .nucleotides import _matching_definitions
    records, status = _matching_definitions(code, set(observed))
    if status == 'unmapped':
        return None, 'unmapped_component'
    if status == 'ambiguous':
        return None, 'ambiguous_contact_atom_mapping'
    analog = records[0].base_analog.upper()
    if analog not in CONTACT_ATOMS:
        return None, 'no_contact_atom_definition'
    positions = dict(records[0].positions)
    names = {k: positions[k] for k in CONTACT_ATOMS[analog] if positions.get(k) not in (None, '', '/')}
    return (names, analog), None


def base_nodes(ctx, node_atoms):
    """C1' atoms for the given O3' node atoms: ([(O3' index, C1' atom)], [missing residue])."""
    found, missing = [], []
    for index, atom in enumerate(node_atoms):
        c1 = ctx.all_atoms.get(atom.residue_key, {}).get(BASE_ATOM, False)
        if c1 is None:
            raise ConversionError(f"Duplicate C1' atom in {atom.description}; choose --altloc")
        if c1 is False:
            missing.append(_residue_identity(atom))
        else:
            found.append((index, c1))
    return found, missing


def _plane_angles(direction, normals):
    length = _norm(direction)
    if length <= 1e-8:
        return None
    return [math.degrees(math.asin(min(1., abs(sum(d*n for d, n in zip(direction, normal)))
                                       / (length*_norm(normal))))) for normal in normals]


def find_base_pairs(ctx, base_node_of, *, hydrogen_cutoff=4.0, max_plane_angle=30.0, min_atoms_per_base=2):
    """Return (stitched edges, candidate records, summary) for C1'-C1' base pairs.

    base_node_of maps residue key -> C1' node index in ctx. Each candidate is one
    residue-pair symmetry orbit with its hydrogen contacts. Significant pairs
    whose residues both have C1' nodes become edges from the first residue's C1'
    (stored position) to the second's in image g.
    """
    hydrogen_cutoff = _positive(hydrogen_cutoff, 'hydrogen_cutoff')
    if max_plane_angle is not None:
        if (isinstance(max_plane_angle, bool) or not isinstance(max_plane_angle, (int, float))
                or not 0 <= max_plane_angle <= 90):
            raise ConversionError('max_plane_angle must be between 0 and 90 degrees, or None')
        max_plane_angle = float(max_plane_angle)
    if isinstance(min_atoms_per_base, bool) or not isinstance(min_atoms_per_base, int) or min_atoms_per_base < 1:
        raise ConversionError('min_atoms_per_base must be a positive integer')
    frame = ctx.frame
    planes, diagnostics = base_planes(ctx, ctx.nucleotide)
    keys = list(planes)
    order = {k: i for i, k in enumerate(ctx.order)}
    contacts, letters = {}, {}
    for key in keys:
        atoms = ctx.all_atoms[key]
        definition, issue = _contact_names(planes[key]['atom'].residue, atoms)
        detail = dict(residue=planes[key]['atom'].residue, description=_residue_identity(planes[key]['atom'])['description'])
        if issue is not None:
            diagnostics.append(dict(detail, issue=issue))
            continue
        names, letters[key] = definition
        missing = sorted(k for k, n in names.items() if n not in atoms)
        if missing or len(names) < len(CONTACT_ATOMS[letters[key]]):
            diagnostics.append(dict(detail, issue='missing_contact_atoms',
                                    atoms=sorted(set(missing) | (set(CONTACT_ATOMS[letters[key]]) - set(names)))))
        chosen = []
        for name in names.values():
            atom = atoms.get(name)
            if name in atoms and atom is None:
                raise ConversionError(f"Base contact atom {name} has alternate locations in "
                                      f"{detail['description']}; choose --altloc")
            if atom is not None and atom.element.upper() in CONTACT_ELEMENTS:
                chosen.append((normalized_name(atom.name), atom, frame.fractional(atom.xyz)))
        contacts[key] = chosen
    covalent = set()
    from .converter import covalent_connections
    for a, b in covalent_connections(ctx.pdb):
        covalent.add(frozenset((a, b)))

    def canonical(a, na, b, nb, g):
        """One description per atom-pair orbit: lowest (residue, atom, operation) first."""
        forward = ((order[a], na, order[b], nb, op_key(g)), (a, na, b, nb, g))
        gi = g.inverse()
        backward = ((order[b], nb, order[a], na, op_key(gi)), (b, nb, a, na, gi))
        return min(forward, backward, key=lambda item: item[0])

    atom_contacts = {}

    def record(a, na, b, nb, g, distance, angles, origin):
        token, (a, na, b, nb, g) = canonical(a, na, b, nb, g)
        if token in atom_contacts:
            atom_contacts[token]['origins'].add(origin)
            return
        atom_contacts[token] = dict(a=a, b=b, g=g, name_a=na, name_b=nb, distance=distance,
                                    angles=angles, origins={origin})

    reach = {k: max((frame.length([x-y for x, y in zip(f, frame.fractional(planes[k]['centroid']))])
                     for _, _, f in contacts[k]), default=0.) for k in contacts}
    hkeys = [k for k in keys if contacts.get(k)]
    for ia, a in enumerate(hkeys):
        fa = frame.fractional(planes[a]['centroid'])
        for b in hkeys[ia:]:
            screen = hydrogen_cutoff + reach[a] + reach[b]
            for op, shift, _ in frame.images(ctx.operations, fa, frame.fractional(planes[b]['centroid']), screen):
                g = _image_op(op, shift)
                identity = g.triplet() == IDENTITY
                if identity and a[0] == b[0]:
                    continue  # same chain copy (NAtopo's interchain-only hydrogen rule)
                _, normal_b = _transform_plane(frame, g, planes[b])
                for na, atom_a, _ in contacts[a]:
                    for nb, atom_b, fb in contacts[b]:
                        image = frame.cartesian(g.apply_to_xyz(fb))
                        distance = math.dist(atom_a.xyz, image)
                        if distance >= hydrogen_cutoff:
                            continue
                        if identity and frozenset((atom_a.serial, atom_b.serial)) in covalent:
                            continue
                        angles = _plane_angles([y-x for x, y in zip(atom_a.xyz, image)],
                                               (planes[a]['normal'], normal_b))
                        if angles is None or (max_plane_angle is not None
                                              and max(angles) > max_plane_angle + 1e-10):
                            continue
                        record(a, na, b, nb, g, distance, angles, 'inferred')

    # Deposited same-image hydrogen annotations (mmCIF conn_type_id hydrog).
    deposited = 0
    if ctx.pdb.conect_types:
        by_serial = defaultdict(list)
        for key in keys:
            for atom in ctx.all_atoms[key].values():
                if atom is not None:
                    by_serial[atom.serial].append((key, atom))
        for (s1, s2), kind in zip(ctx.pdb.conect, ctx.pdb.conect_types):
            if kind.lower() != 'hydrog':
                continue
            for a, atom_a in by_serial.get(s1, ()):
                for b, atom_b in by_serial.get(s2, ()):
                    if a == b:
                        continue
                    angles = _plane_angles([y-x for x, y in zip(atom_a.xyz, atom_b.xyz)],
                                           (planes[a]['normal'], planes[b]['normal']))
                    record(a, normalized_name(atom_a.name), b, normalized_name(atom_b.name), IDENTITY_OP,
                           math.dist(atom_a.xyz, atom_b.xyz), angles, 'file')
                    deposited += 1

    # Group atom contacts per residue-pair orbit.
    grouped = {}
    for c in atom_contacts.values():
        a, b, g, na, nb = c['a'], c['b'], c['g'], c['name_a'], c['name_b']
        gi = g.inverse()
        if (order[b], order[a], op_key(gi)) < (order[a], order[b], op_key(g)):
            a, b, g, na, nb = b, a, gi, nb, na
            c = dict(c, angles=None if c['angles'] is None else c['angles'][::-1])
        token = (order[a], order[b], op_key(g))
        grouped.setdefault(token, dict(a=a, b=b, g=g, contacts=[]))['contacts'].append(dict(c, name_a=na, name_b=nb))

    stitched, records = [], []
    paired = Counter()
    for token in sorted(grouped):
        a, b, g, cs = grouped[token]['a'], grouped[token]['b'], grouped[token]['g'], grouped[token]['contacts']
        left, right = {c['name_a'] for c in cs}, {c['name_b'] for c in cs}
        count = len(cs)
        involution = a == b and g.combine(g).triplet() == IDENTITY
        if involution:
            count = sum(1 if c['name_a'] == c['name_b'] else 2 for c in cs)
            left = right = left | right
        origins = set().union(*(c['origins'] for c in cs))
        origin = 'both' if origins == {'file', 'inferred'} else ('file' if 'file' in origins else 'inferred')
        significant = min(len(left), len(right)) >= min_atoms_per_base
        distances = [c['distance'] for c in cs]
        source, target = base_node_of.get(a), base_node_of.get(b)
        status = 'stitched'
        if not significant:
            status = 'weak_contact'
        elif source is None or target is None:
            status = 'no_graph_node'
        atom_a, atom_b = planes[a]['atom'], planes[b]['atom']
        rec = dict(kind='base_pair', status=status, method=BASE_PAIR_METHOD, significant=significant,
                   origin=origin, atom_pairs=count, atoms_1=sorted(left), atoms_2=sorted(right),
                   base_1=letters.get(a), base_2=letters.get(b),
                   distance_min_A=min(distances), distance_max_A=max(distances),
                   operation=_wrapped(g).triplet(), lattice_translation=_shift(g), image_operation=g.triplet(),
                   residue1=_residue_identity(atom_a), residue2=_residue_identity(atom_b),
                   source=ctx.node_labels[source] if source is not None else None,
                   target=ctx.node_labels[target] if target is not None else None,
                   contacts=[dict(atom_1=c['name_a'], atom_2=c['name_b'], distance_A=c['distance'],
                                  plane_angles_deg=c['angles'], origin=sorted(c['origins']))
                             for c in sorted(cs, key=lambda c: (c['distance'], c['name_a'], c['name_b']))])
        if status == 'stitched':
            c1_a, c1_b = ctx.node_atoms[source], ctx.node_atoms[target]
            image = frame.cartesian(g.apply_to_xyz(frame.fractional(c1_b.xyz)))
            evidence = {k: rec[k] for k in ('origin', 'atom_pairs', 'atoms_1', 'atoms_2', 'base_1', 'base_2',
                                            'distance_min_A', 'distance_max_A')}
            edge = ctx.add_edge(source, target, g, kind='base_pair', distance=math.dist(c1_a.xyz, image),
                                source_atom=c1_a, target_atom=c1_b, evidence=evidence, residues=(a, b))
            if edge is None:
                rec['status'] = 'duplicate_edge'
            else:
                stitched.append(edge)
                paired[a] += 1
                if b != a:
                    paired[b] += 1
        records.append(rec)
    summary = dict(method=BASE_PAIR_METHOD,
                   source='Ported from NAtopo nucleotide_interactions.find_hydrogen_bond_candidates and '
                          'structural_features.FeatureGraph significant hydrogen pairs',
                   base_atom=BASE_ATOM, hydrogen_cutoff_A=hydrogen_cutoff, max_plane_angle_deg=max_plane_angle,
                   min_atoms_per_base=min_atoms_per_base,
                   contact_atoms={k: list(v) for k, v in CONTACT_ATOMS.items()},
                   contact_atom_source='NAtopo bases.xlsx Atom ID columns (Donor/Acceptor ignored)',
                   interchain_only=True, deposited_hydrogen_links=deposited,
                   image_search='every symmetry operation and every lattice translation within reach',
                   terminal_stack_nodes="O3'",
                   base_planes=len(planes), plane_diagnostics=diagnostics,
                   candidate_status_counts=dict(sorted(Counter(r['status'] for r in records).items())),
                   multiple_partner_residues=sum(1 for n in paired.values() if n > 1))
    return stitched, records, summary
