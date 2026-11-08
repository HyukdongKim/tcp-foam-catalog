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

See `data/tcp_foam_catalog.csv` and `figures/tcp_cost_vs_faces.png`. The first ten distinct topologies:

| # | Structure | Source | Space group | Cells | Z12 : Z14 : Z15 : Z16 | c∞ |
|---|---|---|---|---|---|---|
| 1 | A15 (Weaire–Phelan) | 2010 Phys27 | Pm-3n | 8 | 2 : 6 : 0 : 0 | 5.28835 |
| 2 | **N56** | 2010 | I4₁/amd | 76 | 24 : 44 : 8 : 0 | 5.29403 |
| 3 | σ | 2010 Phys26 | P4₂/mnm | 30 | 10 : 16 : 4 : 0 | 5.29421 |
| 4 | Phys25 | 2010 | Cmmm | 30 | 10 : 16 : 4 : 0 | 5.29550 |
| 5 | E11 | 2010 | I4₁/amd | 44 | 16 : 20 : 8 : 0 | 5.29821 |
| 6 | ν (MnSi) | alloy | Immm | 186 | 74 : 80 : 20 : 12 | 5.30141 |
| 7 | N53 | 2010 | P4₂/mmc | 36 | 14 : 14 : 8 : 0 | 5.30152 |
| 8 | N54 | 2010 | P4₂/mmc | 36 | 14 : 16 : 4 : 2 | 5.30184 |
| 9 | E09 | 2010 | Cmc2₁ | 56 | 24 : 20 : 8 : 4 | 5.30479 |
| 10 | Struct21 | 2012 | Cmmm | 42 | 18 : 16 : 4 : 4 | 5.30518 |

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

1. All values are for fixed topology, with no T1 moves. Structures above about 5.40 fill space correctly but are
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
