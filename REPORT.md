# Report: equal-volume foam costs of enumerated TCP structures, and N56

## 1. Problem and quantity

We study periodic partitions of space into cells of equal volume, minimising the total film area (the Kelvin
problem). The figure of merit is the dimensionless cost per cell,

    cost = 2A / (N · V^(2/3)),   V = 1,

where A is the film area in one unit cell (each film counted once) and N is the number of cells in that unit cell.
Reference values are Kelvin 5.3056, σ 5.29423 and Weaire–Phelan (A15) 5.28834. The last two are A from Reddy et al.
(2018) multiplied by (36π)^(1/3).

## 2. Structures

- **84** TCP "space fullerenes" with ≤ 20 cells per fundamental domain (Dutour Sikirić, Delgado-Friedrichs, Deza 2010),
  file `ListDATApaper/*.cgd` in Zenodo 13691342.
- **37** with 21 cells per fundamental domain (Dutour Sikirić, Deza 2012), file `Data4struct/CGDfileBis21/*.cgd`.
- **9** alloy structure files (`SpecialStructure/AlloysStructures`). Five of them repeat a topology already in the 84
  (A15, σ, δ, μ, K₇Cs₆). That leaves **125 distinct topologies**.

Each `.cgd` file gives the cell, space group, the asymmetric cell centres with their coordination (12, 14, 15, 16)
and the neighbour list. We expand these by symmetry with gemmi and keep the explicit neighbour list. For 9 of the
84 structures, the given coordinates do not reproduce the topology as an equal-volume Voronoi/power diagram, so the
neighbour list is necessary. The result is `data/tcp_all.json`: fractional sites, cell, coordination and neighbour
sets.

## 3. Methods

**Flat foam from the graph** (`code/tcpbuild.py`). In a TCP structure every foam vertex is a tetrahedron of four
mutually adjacent cells. The face between cells i and j is the ring of their common neighbours. Each vertex is
placed at the centroid of its tetrahedron. The topology is therefore exactly the enumerated one.

**Own relaxer** (`code/relax.py`). Each face is fan-triangulated about its centre and subdivided r times per edge
(triangles ∝ r²). All mesh points move freely. The objective is the total area. Every cell volume is held at 1
through an augmented Lagrangian, minimised with L-BFGS-B, and the 3 × 3 lattice is relaxed as well
(`relax_box`). Merging of periodic copies uses a 1e-5 tolerance, because the `.cgd` coordinates carry
rounding noise of about 1e-6.

**Continuum limit.** With r = 2, 3, 4 we fit c(r) = c∞ + a/r² + b/r⁴. With only r = 2, 3 we use
c∞ = (9c₃ − 4c₂)/5.

**Surface Evolver 2.70** (`surface_evolver/`). The setup is TORUS_FILLED, every body at volume 1, and periods fixed
at the lattice found by our relaxer. The protocol is in `run.cmd`. At each level it runs `g 100; u; V; g 100`,
then up to 80 rounds of `g 100; u; hessian_seek`, then refines (`r`), for levels 0–3. We extrapolate with
c∞ = (4c₃ − c₂)/3. Surface Evolver itself drops one redundant volume constraint in TORUS_FILLED mode. Removing a
second one by hand left a body unconstrained, so all bodies keep `volume 1`.

## 4. Results

### 4.1 Headline

| Structure | Own relaxer r = 1 / 2 / 3 / 4 | c∞ (own) | Surface Evolver levels 0 / 1 / 2 / 3 | c∞ (SE) |
|---|---|---|---|---|
| WP (A15) | 5.29043 / 5.28881 / 5.28854 / 5.28846 | 5.288345 | 5.29047 / 5.28886 / 5.28846 / 5.28837 | 5.288339 |
| N56 | 5.29595 / 5.29447 / 5.29420 / 5.29412 | 5.294032 | 5.29595 / 5.29447 / 5.29412 / 5.29403 | 5.294000 |
| σ | 5.29606 / 5.29464 / 5.29438 / 5.29430 | 5.294214 | 5.29609 / 5.29465 / 5.29431 / 5.29422 | 5.294193 |

The two programs differ by at most 3e-5 after extrapolation. The N56 − σ gap is −0.00018 (own) and −0.00019 (SE).

### 4.2 Catalogue

See `data/tcp_foam_catalog.csv` and `figures/tcp_cost_vs_faces.png` (the figure shows the 125 v1.0 topologies). The first ten distinct topologies (v1.1, 128 topologies):

| # | Structure | Source | Space group | Cells | Z12 : Z14 : Z15 : Z16 | c∞ |
|---|---|---|---|---|---|---|
| 1 | A15 (Weaire–Phelan) | 2010 Phys27 | Pm-3n | 8 | 2 : 6 : 0 : 0 | 5.28835 |
| 2 | **N56** | 2010 | I4₁/amd | 76 | 24 : 44 : 8 : 0 | 5.29403 |
| 3 | σ | 2010 Phys26 | P4₂/mnm | 30 | 10 : 16 : 4 : 0 | 5.29421 |
| 4 | K (Shoemaker), v1.1 | dataset `SpecialStructure` | Pmmm | 82 | 28 : 42 : 12 : 0 | 5.29525 |
| 5 | Phys25 | 2010 | Cmmm | 30 | 10 : 16 : 4 : 0 | 5.29550 |
| 6 | F (Shoemaker), v1.1 | dataset `SpecialStructure` | P6/mmm | 52 | 18 : 26 : 8 : 0 | 5.29584 |
| 7 | E11 | 2010 | I4₁/amd | 44 | 16 : 20 : 8 : 0 | 5.29821 |
| 8 | J (Shoemaker), v1.1 | dataset `SpecialStructure` | Pmmm | 22 | 8 : 10 : 4 : 0 | 5.29904 |
| 9 | ν (MnSi) | alloy | Immm | 186 | 74 : 80 : 20 : 12 | 5.30141 |
| 10 | N53 | 2010 | P4₂/mmc | 36 | 14 : 14 : 8 : 0 | 5.30152 |

Extrapolation uncertainty is about ±5e-5 for the 22 structures run to r = 4, and up to ±1.3e-4 for the others
(the largest difference between the r2–3 and r2–4 estimates). The same topology given in two different files
agrees within 2e-5 (A15 and Cr₃Si: 5.288347 and 5.288332).

## 5. Verification (WP, σ, N56)

| Check | How | Result |
|---|---|---|
| Convergence | own relaxer r = 1–4 | N56 − σ: −0.00011, −0.00017, −0.00017, −0.00018 |
| Independent program | Surface Evolver levels 0–3 | N56 − σ: −0.00014, −0.00018, −0.00019, −0.00019 |
| Independent start | power-diagram start against graph start (r = 2) | agree within 2e-6 |
| Stability | perturb every mesh point (σ = 0.02) and the lattice (1 %), re-relax, 2 seeds | returns within 2e-6 of the original value |
| Plateau's laws | film angles along triple lines at r = 3 (7,848 for N56) | mean 120.0°, RMS 1.3°, same as WP and σ, falling with r |
| Volumes | all cells | error < 1e-6 (own), < 1e-8 (SE; N56 and WP < 1e-13) |
| Space filling | winding number at 1,500 random points (WP, N55, N57, N58, N60, Struct37) | every point inside exactly one cell |

## 6. Prior work checked

Equal-volume foam areas have been published for A15, σ, H, Z, C14, C15, P, δ, pσ, M, μ, Bergman, S-M and
quasicrystal approximants (Weaire & Phelan 1994; Phelan, Weaire & Brakke 1995; Phelan 1996; Kusner & Sullivan 1996;
Cox et al. 2017; Reddy et al. 2018), plus non-TCP counter-examples (Gabbrielli 2009). The enumeration papers
(Dutour Sikirić et al. 2010, 2012) do not compute areas. We found no published area for N56 or for most of the
enumerated structures. We could not check unpublished data.

## 7. Limitations

1. All catalogue values are for fixed topology. Section 9.2 tests T1 neighbours of WP, σ and N56. Structures above about 5.40 fill space correctly but are
   probably T1-unstable.
2. In Surface Evolver the periods were fixed, not optimised. A relaxed box can only lower the value.
3. Every Surface Evolver level stopped at the 80-round cap, so its values are slight upper bounds. N56 has more
   cells and is likely the least converged, which biases against N56.
4. The σ and N56 Surface Evolver runs started from the graph-built flat foam. The WP run started from a power
   diagram.
5. Struct19 and Struct20 (mean faces 13.81) were borderline at r = 2 (flipped-area fraction about 1e-6, cost about
   6.06). Only their r = 2 values are listed.
6. The mean-faces-above-13.5 structures (N57–N60, Struct19/20/37) relax to costs of 5.41–6.06.

## 8. Reproduction

    pip install -r requirements.txt
    python code/rigor.py              # convergence to r = 4, power-diagram start, perturbation, Plateau angles
    python code/catalog.py 10         # whole catalogue, about 47 core-hours
    python code/space_filling_check.py n20:N56
    evolver -f surface_evolver/run.cmd surface_evolver/n56.fe

Archived outputs: `results/` (`tcp_r1`, `tcp_r2`, `tcp_verify` for r = 3, `rigor`, `catalog`) and
`surface_evolver/*.levels.txt`. Hardware: Apple M4 (10 cores); Python 3.12, numpy 2.5, scipy 1.18.

## 9. Additions in v1.1

### 9.1 Three more TCP structures: Shoemaker F, J, K

The Dutour Sikirić dataset also contains `SpecialStructure/F-phase.cgd`, `J-phase.cgd` and `K-phase.cgd`
(pseudo-hexagonal phases described by Shoemaker & Shoemaker). They are built with `code/shoemaker_build.py`
(`data/tcp_shoemaker.json`) and relaxed at r = 2, 3 with the catalogue settings (`code/shoemaker_run.py`,
archived output `results/shoemaker_phases.jsonl`).

| Structure | Cells | mean faces | r = 2 | r = 3 | c∞ | above hull (9.3) |
|---|---|---|---|---|---|---|
| K | 82 | 13.4634 | 5.295690 | 5.295443 | 5.29525 | +0.00031 |
| F | 52 | 13.4615 | 5.296294 | 5.296044 | 5.29584 | +0.00048 |
| J | 22 | 13.4545 | 5.299477 | 5.299233 | 5.29904 | +0.00210 |

K is now fourth, below Phys25. With them, 13 of 128 distinct topologies beat Kelvin's foam. Four other alloy
files in `SpecialStructure/AlloysStructures` (R, X, a second ν file, H-complex) could not be used: their cells or
coordinates do not reproduce a TCP net (for example the R-phase cell has angles 74.57°, 90°, 90° and still gives
63 instead of 53 sites after correction).

### 9.2 Non-TCP topologies by T1 moves

A 3D foam with Plateau vertices is dual to a periodic triangulation: cells ↔ vertices, faces ↔ edges, foam edges ↔
triangles, foam vertices ↔ tetrahedra. A T1 move is a Pachner 2–3 flip (a foam edge becomes a triangular face) or its
inverse 3–2. `code/trifoam.py` stores the triangulation (tetrahedra up to lattice translation), applies flips, checks
validity (each triangle in two tetrahedra, every edge link one cycle, every vertex link a sphere) and builds faces
for the relaxer. Rebuilding WP and σ from the triangulation reproduces the graph-built values (r = 1: 5.290432 for
WP, identical to 1e-8). These topologies include 3-, 4- and 7-gon faces and need not be power diagrams.

**Beam search** (`code/t1search.py`, r = 1, box relaxed, results `results/t1/t1_*.jsonl.gz`): 9,808 topologies
from WP (1×1×1 to depth 8, 2×1×1 to depth 6, 2×2×1 to depth 3), σ (depth 3) and N56 (depth 1). None was lower than
its start. Best Δ (r = 1): WP 1×1×1 +0.0079, 2×1×1 +0.0041, 2×2×1 +0.0022, σ +0.0018, N56 +0.0009, i.e. a
dilute single T1 costs about 0.054–0.069 (Δ × N). Cost grows about linearly with the number of flips.

**Caution.** These best candidates all carry one triangular face, and on relaxation that face shrinks to zero
area (`code/t1_diag.py`): the foam is trying to undo the T1. Their Δ therefore falls with resolution towards the
parent value (WP 2×2×1: +0.00215, +0.00167, +0.00012 at r = 1, 2, 3), and the r2–3 extrapolation of §3 is invalid
for them; applied blindly it gives a spurious value below WP. They are not new stable foams.

**Triangle-free topologies** (`code/t1_stable.py`, results `results/t1/t1s_*.jsonl.gz`): exhaustive BFS of flip
sequences, deduplicated by a topological fingerprint, without relaxing intermediates; only end states without
triangular faces are relaxed. (v1.1.0 said such states occur only at odd depth; that is wrong: WP 1×1×1 has
6 triangle-free states at depth 3 and 6 more at depth 4.)

| Start | depth | triangle-free topologies | best Δ (r = 1) | Δ × N |
|---|---|---|---|---|
| WP 1×1×1 | 4 | 12 | +0.0264 | 0.21 |
| WP 2×1×1 | 4 (300,000-state cap) | 8 | +0.0137 | 0.22 |
| σ 1×1×1 | 3 | 37 | +0.0065 | 0.19 |

The smallest stable non-TCP defect (flips 2–3, 2–3, 3–2; five 4-gons) costs about 0.2, three times a collapsed T1.
WP, σ and N56 are local minima among topologies within three to four T1 moves. Limits: the fingerprint is not a
full isomorphism test (distinct topologies may be merged), the number of cells is fixed (no T2), depth 5 was
not reached, and the WP 2×1×1 depth-4 level was cut by the state cap, so its depth-4 triangle-free states are
incomplete.

### 9.3 Lower convex hull of cost against mean faces

Plotting c∞ against mean faces per cell Z̄ for all 128 topologies, the lower convex hull has only three vertices:
Z̄ = 13.333 (C15 family, 5.32421), 13.467 (σ, 5.29421) and 13.5 (A15, 5.28833). Every other structure lies above
it (N56 +0.0011, K +0.0003). Structures with Z̄ between σ and A15 are A15/σ intergrowths such as N56, and their
position above the hull indicates a positive interface cost. A larger TCP would have to fall well below this hull to
beat WP. This is an empirical trend, not a proof.
