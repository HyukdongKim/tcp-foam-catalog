# Checking and challenging these results

Challenges are welcome. The aim of this repository is to let anyone confirm or overturn the numbers.

## Most useful reports

1. **Prior computation.** You know of a published or unpublished value for N56, or for another enumerated
   structure in `data/tcp_foam_catalog.csv`. Please give the reference.
2. **Numerical disagreement.** Your own computation (Surface Evolver, another code, or ours with other settings)
   gives a different value or a different order. Please state the structure, the value, the mesh or resolution, how
   you extrapolated, and whether the unit cell was relaxed.
3. **Topology or data error.** A structure in `data/tcp_all.json` does not match its source `.cgd`, a name is
   wrong, or two entries are the same topology.
4. **Method concerns.** Fixed topology, fixed periods in Surface Evolver, the extrapolation model, the validity
   criteria, or anything else you think could change a conclusion.

Please use the issue templates. If your report changes a number, we will record it in the changelog and keep the
original value visible.

## Running the checks

    pip install -r requirements.txt
    python code/rigor.py                       # writes rerun/rigor.jsonl
    python code/catalog.py 10                  # writes rerun/catalog.jsonl (resumable)
    python code/space_filling_check.py n20:N56
    python code/export_fe.py data/tcp_all.json n20:N56 - n56.fe --voro --fixv   # a fresh Surface Evolver file

`export_fe.py` writes the flat foam in the cell stored in `tcp_all.json`, or in a relaxed cell given as JSON
(`{"n20:N56": [[...],[...],[...]]}`) in place of `-`. The archived Surface Evolver inputs in `surface_evolver/` are the
exact files that produced the published levels.
