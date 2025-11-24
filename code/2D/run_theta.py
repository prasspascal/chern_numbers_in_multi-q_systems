import json
import os
import sys
import time

import numpy as np
from constants import deg
from input import params_thetas_list as params_list
from model import construct_hilbert_space, spectrum
from mpicore import MPIControl
from qspace import QSpace

# ----------------------------------------------------------------------
# setup mpi
# ----------------------------------------------------------------------

mpiv = MPIControl()

for params in params_list:
    mpiv.start_clock()

    # ----------------------------------------------------------------------
    # generate output outdirectory
    # ----------------------------------------------------------------------

    outdir = ""
    if mpiv.is_root():

        # unique time stamp
        time.sleep(1)
        stamp = str(round(time.time()))

        # output outdirectory
        outdir = "./out/" + stamp
        os.mkdir(outdir)

        with open(outdir + "/params.json", "w") as f:
            json.dump(params, f)

    outdir = mpiv.comm.bcast(outdir, root=0)

    # ----------------------------------------------------------------------
    # calculate the spectrum
    # ----------------------------------------------------------------------

    # determine admissible q-vectors
    Q = QSpace(params["system_sizes"])

    qs, ns = Q.get_qlist()
    n_q = len(qs)

    if mpiv.is_root():

        np.save(outdir + "/qs.npy", qs)
        np.save(outdir + "/ns.npy", ns)

    for i in range(n_q):
        if mpiv.my_turn(i):
            print("current index:", i + 1, n_q)
            sys.stdout.flush()

            labels, states = construct_hilbert_space(ns[i], ns[i])

            eigvals = 2 * (ns[i] * ns[i])

            spec = np.zeros(eigvals, dtype=np.float64)
            spec = spectrum(
                states,
                ns[i],
                ns[i],
                qs[i],
                0.0,
                np.array([0,0]),
                params["shift"],
                params["mag"],
                params["t"],
                params["m"],
                params["texture"],
                "periodic",
            )

            key = str(i).zfill(4)
            np.save(outdir + "/spec_" + key + ".npy", spec)

    # ---------------------------------------------------------------------
    # finalize
    # ----------------------------------------------------------------------

    mpiv.barrier()
    mpiv.print("Done!")
    mpiv.stop_clock()

    if mpiv.is_root():

        walltime = mpiv.get_time()
        cputime = mpiv.size * mpiv.get_time() / 3600.0

        mpiv.print("Walltime: ", walltime)
        mpiv.print("CPUs: ", mpiv.size * mpiv.get_time())
        mpiv.print("outdir: ", outdir)
        np.savetxt(outdir + "/cputime.txt", [cputime])

mpiv.finalize()