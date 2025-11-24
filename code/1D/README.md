# Python code

The Python scripts `./run_*` generate the data used in the thesis from the input parameters stored in `./input.py`.
The output is sent to the directory `./out`.
The parameters for each generated data set are recorded in the respective `./out/*/params.json` files.

The Python files `./model.py` and `./constants.py` contains all methods used to set up the spin Hamiltonian, construct the Fermi projections, and evaluate the Chern numbers.

The Python file `./qspace.py` contains the methods to sample rational pitch parameters.

The Python file `./kpm.py` contains the methods for the Kernel Polynomial Method.
See DOI: [10.1103/RevModPhys.78.275](https://doi.org/10.1103/RevModPhys.78.275) for a detailed description of the method itself.

The Python file `./mpicore.py` sets up MPI for Python to provide Python bindings for the Message Passing Interface (MPI) standard, allowing Python applications to exploit multiple processors on workstations, clusters, and supercomputers.

| Scripts | Figures and Tables |
|---|---|
| `./run_theta.py` | Fig. 5.1 |
| `./run_theta_kpm.py` `./run_chern_fermi.py` | Fig. 5.4 |
| `./run_chern.py` | Tbl. 5.1 |
