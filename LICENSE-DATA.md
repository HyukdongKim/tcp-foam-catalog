# License for data, results, figures and documentation

Everything in this repository outside `code/` (including `data/`, `results/`, `surface_evolver/`, `figures/`,
`README*.md`, `REPORT.md`) is licensed under the
**Creative Commons Attribution 4.0 International License (CC BY 4.0)**:
https://creativecommons.org/licenses/by/4.0/

Copyright (c) 2026 Hyukdong Kim.

## Third-party material

`data/source_cif/*.cif` are unmodified files from the Crystallography Open Database (https://www.crystallography.net/,
COD 2310299, 2310312, 2311049). The original data were provided by IUCr Journals; please cite the journal articles
named in each file. `data/source_cif/M_phase_NbNiAl_Shoemaker1967.cif` was written by the author from Table 2 of
C. B. Shoemaker & D. P. Shoemaker, Acta Cryst. 23 (1967) 231. `data/tcp_cif.json` is derived from these files.

`data/source_cgd/*.cgd` are unmodified files, and `data/tcp_all.json` and the topologies in
`surface_evolver/*.fe` are adapted, from:

> M. Dutour Sikirić (2024), *Dataset and programs for the work on space-fullerenes*, Zenodo, record 13691342,
> https://zenodo.org/records/13691342, licensed CC BY 4.0.

Changes: symmetry expansion of cell centres and neighbour lists (gemmi), unit-volume scaling, conversion to JSON and
to Surface Evolver datafiles.
