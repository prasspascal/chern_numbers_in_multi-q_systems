import json
import os
import time
import sys

import numpy as np
from constants import pi
from input import params_fermi_list as params_list
from model import construct_hilbert_space, spectrumV, spec_projection, phason_derivative, translation_derivative, trace
from mpicore import MPIControl

# ----------------------------------------------------------------------
# setup mpi
# ----------------------------------------------------------------------

mpiv = MPIControl()

for j,params in enumerate(params_list):
    mpiv.start_clock()

    # size of unit cell and texture parameter
    n = params["system_sizes"]
    q = params["q"]

    fermis = np.linspace(params["min_fermi"],params["max_fermi"],params["n_fermi"])

    tau1u1_partial = np.zeros(0,dtype=np.complex128)
    partial = np.zeros(0)
    
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


    # ----------------------------------------------------------------------
    # calculate Chern numbers
    # ----------------------------------------------------------------------

    # calculate the spectrum
    labels, states = construct_hilbert_space(n)

    eigenvals , eigenvecs = spectrumV(
                            states,
                            n,
                            q,
                            0.0,
                            2*q,
                            0.0,
                            1,
                            params["t"],
                            params["m"],
                            params["texture"],
                            True,
                            )

    for i,fermi in enumerate(fermis):
        if mpiv.my_turn(i):
            print(i,"fermi: ",fermi)
            sys.stdout.flush()

            # list for partial derivatives of the spectral operator
            dP = np.zeros((2,2*n,2*n),dtype=np.complex128)

            # calculate the spectral projection operator
            P = spec_projection(eigenvals,eigenvecs,fermi,params["kbT"])

            dP[0] = translation_derivative(P, states, n)

            dP[1] = phason_derivative(
                                    P,
                                    params["delta"],
                                    states,
                                    n,
                                    q,
                                    0.0,
                                    2*q,
                                    0.0,
                                    1,
                                    params["t"],
                                    params["m"],
                                    params["texture"],
                                    True,
                                    fermi,
                                    params["kbT"]
                                    )

            tau1u1_partial = np.append(tau1u1_partial,[(2*pi*1j) * 2*1j*trace(n, P.dot(dP[0].dot(dP[1]))).imag])

            partial = np.append(partial,i)

    mpiv.barrier()
    mpiv.print("Done!")
    mpiv.stop_clock()

    # gather and flatten
    index = mpiv.gather(partial)

    if mpiv.is_root():
        index = np.array([item for sublist in index for item in sublist])
        #mpiv.print(index)
        #mpiv.print(index[np.argsort(index)])

    # gather lists of lists chern numbers
    tau1u1_list = mpiv.gather(tau1u1_partial)

    if mpiv.is_root():
        np.save(outdir + "/spectrum.npy", eigenvals)
        # flatten list of lists to numpy array, sort numpy array by index, and save numpy array
        tau1u1_list = np.array([item for sublist in tau1u1_list for item in sublist])
        tau1u1_list = tau1u1_list[np.argsort(index)]
        np.save(outdir + "/ch_tau1u1.npy", tau1u1_list)

    if mpiv.is_root():

        walltime = mpiv.get_time()

        mpiv.print("Walltime: ", walltime)
        mpiv.print("outdir: ", outdir)
        np.savetxt(outdir + "/walltime.txt", [walltime])

mpiv.finalize()
