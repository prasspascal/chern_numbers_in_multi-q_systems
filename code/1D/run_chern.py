import json
import os
import sys
import time

import numpy as np
# from constants import deg
from input import params_chern_list as params_list
from model import construct_hilbert_space, spectrumV, spec_projection, phason_derivative, translation_derivative, chern_character
from mpicore import MPIControl

# ----------------------------------------------------------------------
# setup mpi
# ----------------------------------------------------------------------

mpiv = MPIControl()

for j,params in enumerate(params_list):
    mpiv.print("Current Set of Parameters: ", j+1)
    time.sleep(1)

    # size of unit cell and texture parameter
    ns = params["system_sizes"]
    qs = params["q"]
    fermis = params["fermi"]
    n_q = len(qs)

    for i in range(n_q):
        if mpiv.my_turn(i):
            print("current index:", i + 1, n_q)
            sys.stdout.flush()

            # ----------------------------------------------------------------------
            # generate output outdirectory
            # ----------------------------------------------------------------------
            outdir = ""
            # unique time stamp
            time.sleep(2*i)
            stamp = str(round(time.time()))

            # output outdirectory
            outdir = "./out/" + stamp
            os.mkdir(outdir)

            with open(outdir + "/params.json", "w") as f:
                tmp =  {key: value for key, value in params.items()}
                tmp["system_sizes"]=[ns[i]]
                tmp["q"]=qs[i]
                tmp["fermi"]=fermis[i]
                json.dump(tmp, f)

            t0 = time.time()

            # ----------------------------------------------------------------------
            # calculate Chern numbers
            # ----------------------------------------------------------------------

            # list for partial derivatives of the spectral operator
            dP = np.zeros((2,2*ns[i],2*ns[i]),dtype=np.complex128)

            # calculate the spectrum
            labels, states = construct_hilbert_space(ns[i])

            eigenvals , eigenvecs = spectrumV(
                                    states,
                                    ns[i],
                                    qs[i],
                                    0.0,
                                    2*qs[i],
                                    0.0,
                                    1,
                                    params["t"],
                                    params["m"],
                                    params["texture"],
                                    True,
                                    )

            # calculate the spectral projection operator
            P = spec_projection(eigenvals,eigenvecs,fermis[i])

            print("current index:", i + 1, n_q,"Spectral Projection: ",  time.time()-t0, np.linalg.norm(P-np.transpose(np.conj(P))) )
            sys.stdout.flush()

            # Calculate the derivatives of the spectral projection operator

            #dP_tau1
            dP[0] = translation_derivative(P, states, ns[i])

            print("current index:", i + 1, n_q,"Translation derivative: ",  time.time()-t0, np.linalg.norm(dP[0]-np.transpose(np.conj(dP[0]))) )
            sys.stdout.flush()

            #dP_phi1
            dP[1] = phason_derivative(
                            P,
                            params["delta"],
                            states,
                            ns[i],
                            qs[i],
                            0.0,
                            2*qs[i],
                            0.0,
                            1,
                            params["t"],
                            params["m"],
                            params["texture"],
                            True,
                            fermis[i]
                            )

            print("current index:", i + 1, n_q,"Phason derivative: ",  time.time()-t0, np.linalg.norm(dP[1]-np.transpose(np.conj(dP[1]))) )
            sys.stdout.flush()


            # Compute the Chern numbers for all even subsets of the labels which are True in J
            chern_numbers = chern_character(P,dP,ns[i])

            np.save(outdir + "/chern_numbers_tau1u1.npy", dict(chern_numbers)) # numba does not return an ordinary dictionary but a numba.typed.Dict

            walltime = time.time()-t0
            print("current index:", i + 1, n_q,"Walltime: ", walltime)
            sys.stdout.flush()
            np.savetxt(outdir + "/walltime.txt", [walltime])

    # ---------------------------------------------------------------------
    # finalize
    # ----------------------------------------------------------------------

    mpiv.barrier()
    mpiv.print("Done!")
    mpiv.stop_clock()

    if mpiv.is_root():

        walltime = mpiv.get_time()

        mpiv.print("Walltime: ", walltime)
        mpiv.print("outdir: ", outdir)

mpiv.finalize()
