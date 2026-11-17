# 8D93 example outputs

Generated with TopoStich 0.4.0 from PDB entry 8D93 (`8d93.cif`, SHA-256
`9f67eb6a7193d1ff06b68105e4ba8a066f4b02bef55787fc096404a99fa1d7bb`; the input
itself is not included). They were added on request as reference outputs.
They are not packaged, and they are regenerated rather than edited. The
`8d93_stitched_hbond_*` files were generated with TopoTrace 0.5.0 from the same
input (see [below](#base-pairs-8d93_stitched_hbond)).

```sh
topostich-cgd 8d93.cif --stitch-nicks --stitch-terminal-stacks --name 8D93_stiched \
  --out 8D93_stiched.cgd --receipt 8D93_stiched.json \
  --cmm 8D93_stiched.cmm --crystal-cmm 8D93_stiched_symmetry.cmm \
  --asu-graphs 8D93_stiched_asu.json \
  --asu-chain-html 8D93_stiched_chains.html --asu-residue-html 8D93_stiched_residues.html \
  --asu-chain-png 8D93_stiched_chains.png --asu-residue-png 8D93_stiched_residues.png
topostich-cgd 8d93.cif --name 8D93_unstiched \
  --out 8D93_unstiched.cgd --receipt 8D93_unstiched.json \
  --cmm 8D93_unstiched.cmm --crystal-cmm 8D93_unstiched_symmetry.cmm \
  --asu-graphs 8D93_unstiched_asu.json \
  --asu-chain-html 8D93_unstiched_chains.html --asu-residue-html 8D93_unstiched_residues.html \
  --asu-chain-png 8D93_unstiched_chains.png --asu-residue-png 8D93_unstiched_residues.png
```

| File | Content |
|---|---|
| `*.cgd` | Systre CRYSTAL graph, one O3′ node per nucleotide |
| `*.json` | Receipt: node identities, edges, stitched joins and every candidate with its status |
| `*.cmm`, `*_symmetry.cmm` | ChimeraX marker views: asymmetric unit and intact symmetry copies |
| `*_asu.json` | ASU quotient graphs at chain and residue resolution, with periodicity |
| `*_chains.html`, `*_residues.html` | Interactive views of those graphs, open in a browser |
| `*_chains.png`, `*_residues.png` | The same graphs as images (needs `topostich[png]`) |

**Unstitched:** 42 nodes and 38 edges in four separate finite strands (A=12,
B=7, C=14, D=9).

**Stitched:** 42 nodes and 42 edges. The four added edges are:

- two backbone nick joins, both with code `3_(0,0,0)`: A:DA112 O3′ to D:DC113 P at
  2.93 Å, and B:DT114 to B:DA105 at 3.12 Å, which joins B to its 3-fold copies in a closed finite ring;
- two 5′-3′ chain-end stacks, both with code `9_(0,0,0)`: C:DC114 on C:DT101 at
  3.38 Å, and D:DA121 on A:DG101 at 3.91 Å.

Two further end stacks were not stitched because those ends are already
nick-joined. The ASU graph gives A+D a 1-periodic rod along [2/3, 1/3, 1/3],
C a rod along [0, 0, 1], and B finite.

Both files pass `tools/check_stitched_cgd.py` with Systre 19.6.0 (the edge
counts match the receipts).

## Base pairs: `8d93_stitched_hbond`

Generated with TopoTrace 0.5.0 (2026-11-17) from the same `8d93.cif`, with C1′
base nodes and base pairs as well as nick and chain-end stack stitching:

```sh
topotrace-cgd 8d93.cif --base-pairs --stitch-nicks --stitch-terminal-stacks --name 8d93_stitched_hbond \
  --out 8d93_stitched_hbond.cgd --receipt 8d93_stitched_hbond.json \
  --cmm 8d93_stitched_hbond.cmm --crystal-cmm 8d93_stitched_hbond_symmetry.cmm \
  --asu-graphs 8d93_stitched_hbond_asu.json \
  --asu-chain-html 8d93_stitched_hbond_chains.html --asu-residue-html 8d93_stitched_hbond_residues.html \
  --asu-chain-png 8d93_stitched_hbond_chains.png --asu-residue-png 8d93_stitched_hbond_residues.png
```

**84 nodes and 105 edges.** Every nucleotide has both an O3′ node
(`S<seg>_N<k>`) and a C1′ base node (`S<seg>_B<k>`). The edges are:

- 38 backbone trace edges;
- 42 O3′–C1′ links, one per nucleotide;
- the same 2 nick joins and 2 chain-end stacks as `8D93_stiched`, which stay on
  O3′ nodes;
- 21 C1′–C1′ base pairs. These are NAtopo's significant hydrogen pairs: at
  least 2 contact atoms on each base, contacts under 4.0 Å and within 30° of
  both base planes. 18 lie inside the ASU (`1_(0,0,0)`). Three join symmetry
  copies: A:DG101 and A:DA102 pair C at `9_(0,0,0)`, and A:DA112 pairs B:DT114 at
  `3_(0,0,0)`. 18 pairs are also deposited `hydrog` annotations
  (origin `both`), and 3 were found only geometrically.

Every base is paired, so there are no unpaired C1′ nodes. One weak contact,
A:DC104 N4…C:DT112 O4 with a single atom on each base, is reported in
`base_pair_candidates` but not linked. These 21 pairs match NAtopo's own
`base_pairs.csv` for 8D93 exactly.

The base pairs tie all four strands into a single **3-periodic framework**
along [1/3, 2/3, 2/3], [0, 1, 0] and [0, 0, 1], at both chain and residue
resolution. Systre 19.6.0 agrees: one 3-periodic component of 252 nodes in the
primitive cell. `tools/check_stitched_cgd.py` passes, with 315 Systre edges
equal to the receipt. Without the base pairs, `8D93_stiched` above has only
1-periodic rods and a finite ring. Teal edges in the views and images are base
pairs. In the CMM files, C1′ nodes are the small markers inside each strand's
model, and the base pairs have their own **base pairs** model.
