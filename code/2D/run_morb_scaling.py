import json
import os
import time
import sys

import numpy as np
from input import params_morb_scaling_list as params_list
from model import construct_hilbert_space, spectrumV, spec_projection, translation_derivative, orbital_magnetisation
from mpicore import MPIControl

# ----------------------------------------------------------------------
# setup mpi
# ----------------------------------------------------------------------

mpiv = MPIControl()

for params in params_list:
    mpiv.start_clock()

    # size of unit cell and texture, field, and system parameter
    n = params["system_sizes"]

    qs = params["q"]
    shift = params["shift"]
    mag   = params["mag"]

    fluxes = params["flux"]

    t = params["t"]
    m = params["m"]

    # gap labels for ids computation
    gap_pos = params["gap_pos"]
    block = params["block"]
    tau1tau2 = params["tau1tau2"]
    tau1tau2u1u2 = params["tau1tau2u1u2"]

    fermi_partial      = np.zeros(0,dtype=np.complex128)
    pre_fermi_partial  = np.zeros(0,dtype=np.complex128)
    post_fermi_partial = np.zeros(0,dtype=np.complex128)

    tau1tau2_partial = np.zeros(0,dtype=np.complex128)
    morb_partial     = np.zeros(0,dtype=np.complex128)
    morbLC_partial   = np.zeros(0,dtype=np.complex128)
    morbIC_partial   = np.zeros(0,dtype=np.complex128)

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
            tmp =  {key: value for key, value in params.items()}
            json.dump(tmp, f)

    outdir = mpiv.comm.bcast(outdir, root=0)
    
    # ----------------------------------------------------------------------
    # calculate Chern numbers
    # ----------------------------------------------------------------------

    # calculate the spectrum
    labels, states = construct_hilbert_space(n, n)

    for i in range(len(qs)*len(fluxes)):
        if mpiv.my_turn(i):

            j = i % len(fluxes)
            k = i //len(fluxes)

            # print(i," theta: ",qs[k],"flux: ",fluxes[j])
            # sys.stdout.flush()

            H, eigenvals , eigenvecs = spectrumV(
                                                states,
                                                n,
                                                n,
                                                qs[k],
                                                0.0,
                                                np.array([0.0,0.0]),
                                                shift,
                                                mag,
                                                t,
                                                m,
                                                params["texture"],
                                                "periodic",
                                                params["eta"],
                                                fluxes[j]
                                                )
            
            # compute ids from gap labels
            ids = block + tau1tau2 * fluxes[j] - tau1tau2u1u2 * qs[k]**2 + 0.5/n**2
            
            # catch ids boundaries where Chern numbers change
            if ids<=0:
                ids = -ids
            elif ids>= 2:
                ids = 4-ids

            print(i," theta: ",qs[k],"flux: ",fluxes[j],"ids: ",ids,"fermi: ",eigenvals[int(ids*n*n)-2:int(ids*n*n)+1])
            sys.stdout.flush()

            # calculate the spectral projection operator
            P = spec_projection(n,n,eigenvecs,ids)

            # list for partial derivatives of the spectral operator
            dP_J = np.zeros((2,2*n*n,2*n*n),dtype=np.complex128)

            dP_J[0] = translation_derivative(P, 0, states, n, n)

            dP_J[1] = translation_derivative(P, 1, states, n, n)

            # fermi energy to given IDS (at gap_pos of the corresponding gap)
            if gap_pos == "top":
                fermi = eigenvals[int(ids*n*n)]
            elif gap_pos == "middle":
                fermi = (eigenvals[int(ids*n*n)] + eigenvals[int(ids*n*n)-1])/2
            else:
                fermi = eigenvals[int(ids*n*n)-1]

            fermi_partial = np.append(fermi_partial,[eigenvals[int(ids*n*n)-1]])
            pre_fermi_partial  = np.append(pre_fermi_partial,[eigenvals[int(ids*n*n)-2]])
            post_fermi_partial = np.append(post_fermi_partial,[eigenvals[int(ids*n*n)]])

            # contributions to orbital magnetisation
            total_morb = orbital_magnetisation(H,P,dP_J,fermi,n,n)

            morb_partial   = np.append(morb_partial,[total_morb[0]])
            morbLC_partial = np.append(morbLC_partial,[total_morb[1]])
            morbIC_partial = np.append(morbIC_partial,[total_morb[2]])

            tau1tau2_partial = np.append(tau1tau2_partial,[total_morb[3]])

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
    fermi_list = mpiv.gather(fermi_partial)
    pre_fermi_list = mpiv.gather(pre_fermi_partial)
    post_fermi_list = mpiv.gather(post_fermi_partial)

    tau1tau2_list = mpiv.gather(tau1tau2_partial)
    morb_list = mpiv.gather(morb_partial)
    morbLC_list = mpiv.gather(morbLC_partial)
    morbIC_list = mpiv.gather(morbIC_partial)

    if mpiv.is_root():
        # flatten list of lists to numpy array, sort numpy array by index, and save numpy array
        fermi_list = np.array([item for sublist in fermi_list for item in sublist])
        fermi_list = fermi_list[np.argsort(index)]
        np.save(outdir + "/fermi.npy", fermi_list)

        pre_fermi_list = np.array([item for sublist in pre_fermi_list for item in sublist])
        pre_fermi_list = pre_fermi_list[np.argsort(index)]
        np.save(outdir + "/pre_fermi.npy", pre_fermi_list)

        post_fermi_list = np.array([item for sublist in post_fermi_list for item in sublist])
        post_fermi_list = post_fermi_list[np.argsort(index)]
        np.save(outdir + "/post_fermi.npy", post_fermi_list)

        tau1tau2_list = np.array([item for sublist in tau1tau2_list for item in sublist])
        tau1tau2_list = tau1tau2_list[np.argsort(index)]
        np.save(outdir + "/ch_tau1tau2.npy", tau1tau2_list)

        morb_list = np.array([item for sublist in morb_list for item in sublist])
        morb_list = morb_list[np.argsort(index)]
        np.save(outdir + "/morb.npy", morb_list)

        morbLC_list = np.array([item for sublist in morbLC_list for item in sublist])
        morbLC_list = morbLC_list[np.argsort(index)]
        np.save(outdir + "/morbLC.npy", morbLC_list)

        morbIC_list = np.array([item for sublist in morbIC_list for item in sublist])
        morbIC_list = morbIC_list[np.argsort(index)]
        np.save(outdir + "/morbIC.npy", morbIC_list)

    if mpiv.is_root():

        walltime = mpiv.get_time()

        mpiv.print("Walltime: ", walltime)
        mpiv.print("outdir: ", outdir)
        np.savetxt(outdir + "/walltime.txt", [walltime])

mpiv.finalize()
