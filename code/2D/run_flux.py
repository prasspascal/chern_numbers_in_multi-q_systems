#import datetime
import json
import os
import sys
import time
#from fractions import Fraction

import numpy as np
from input import params_fluxes_list as params_list
from model import construct_hilbert_space, spectrum
from mpicore import MPIControl
#from numba import jit
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

    if params["bc"] == "periodic":

        # determine admissible q-vectors
        F = QSpace(params["system_sizes"])

        fs, ns = F.get_qlist()
        n_f = len(fs)

        if mpiv.is_root():

            np.save(outdir + "/qs.npy", fs)
            np.save(outdir + "/ns.npy", ns)

        for i in range(n_f):
            if mpiv.my_turn(i):
                print("current index:", i + 1, n_f)
                sys.stdout.flush()

                labels, states = construct_hilbert_space(ns[i], ns[i])

                eigvals = 2 * (ns[i] * ns[i])

                spec = np.zeros(eigvals, dtype=np.float64)
                spec = spectrum(
                    states,
                    ns[i],
                    ns[i],
                    params["q"],
                    0.0,
                    np.array([0,0]),
                    params["shift"],
                    params["mag"],
                    params["t"],
                    params["m"],
                    params["texture"],
                    "periodic",
                    flux=fs[i]
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
