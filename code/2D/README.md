# Python code

The Python scripts `./run_*` generate the data used in our work from the input parameters stored in `./input.py`. The output is sent to the directory `./out`.
The parameters for each generated data set are recorded in the respective `*/params.json` files.

The Python files `./model.py`, `./lattice.py`, and `./constants.py` contains all methods used to set up the spin Hamiltonian, construct the Fermi projections, and evaluate the Chern numbers.

The Python file `./qspace.py` contains the methods to sample rational pitch parameter.

The Python file `./kpm.py` contains the methods for the Kernel Polynomial Method. See DOI: [10.1103/RevModPhys.78.275](https://doi.org/10.1103/RevModPhys.78.275) for a detailed description of the method itself.

The Python file `./mpicore.py` sets up MPI for Python to provide Python bindings for the Message Passing Interface (MPI) standard, allowing Python applications to exploit multiple processors on workstations, clusters, and supercomputers.

| Scripts | Figures and Tables |
|---|---|
| `./run_chern_ids.py`  | Fig. 4.1. |
| `./run_theta.py` | Fig. 5.2 |
| `./run_theta_kpm.py` `./run_chern_fermi.py` | Fig. 5.5 5.6 5.10 5.14 5.15 5.16 5.17 6.1 |
| `./run_chern.py` | Tbl. 5.2 5.3 5.4 |
| `./run_shifts_kpm.py` `./run_mags_kpm.py` | Fig. 5.11 5.12 5.23 5.24 |
| `./run_shifts.py` `./run_mags.py` `./run_chern_shifts_ids.py` `./run_chern_mags_ids.py` | Fig. 5.19 5.20 5.25 5.26 6.2 |
| `./run_chern_diagram.py` | Fig. 5.21 5.22 6.3 |
| `./run_flux.py` | Fig. 6.4 |
| `./run_flux_kpm.py` `./run_chern_fermi.py` | Fig. 6.5 |
| `./run_theta_kpm.py` `./run_flux_kpm.py` | Fig. 6.6 6.10 |
| `./run_theta_kpm.py` `./run_flux_kpm.py` `./run_chern_fermi.py`| Fig. 6.7 6.11 |
| `./run_morb_scaling.py` | Fig. 6.8 |