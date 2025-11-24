from constants import pi
import numpy as np

# ----------------------------------------------------------------------
# setup model
# ----------------------------------------------------------------------

params_thetas_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        "bc": "periodic",  #'open' or 'periodic'
        # -- settings of periodic boundary conditions
        "system_sizes": [i for i in range(19, 80)],
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        "bc": "periodic",  #'open' or 'periodic'
        # -- settings of periodic boundary conditions
        "system_sizes": [i for i in range(19, 80)],
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    # vortex lattice
    {
        "texture": "sdw",
        "bc": "periodic",  #'open' or 'periodic'
        # -- settings of periodic boundary conditions
        "system_sizes": [i for i in range(19, 80)],
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
]

params_thetas_kpm_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.7,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": -0.5,
    },
    # vortex lattice
    {
        "texture": "sdw",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 0/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    # helical skyrmion lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 37/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 71/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 107/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 151/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    # meron-antimeron lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "flux": 71/400, # fixed flux value for Chern number calculation
        # -- settings kpm
        "scale": 12,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": pi/2,
        "mag": 0.0,
    },
]

params_fluxes_list = [
    {
        "texture": "skx",
        "bc": "periodic",  #'open' or 'periodic'
        # -- settings of periodic boundary conditions
        "system_sizes": [i for i in range(19, 80)],
        "q": 0.0,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 0,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
        {
        "texture": "skx",
        "bc": "periodic",  #'open' or 'periodic'
        # -- settings of periodic boundary conditions
        "system_sizes": [i for i in range(19, 80)],
        "q": 0.0,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
]

params_fluxes_kpm_list = [
    # skyrmion lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 0/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 37/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 71/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 107/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    {   
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 151/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
    },
    # meron-antimeron lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 400,
        "q": 71/400,
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": pi/2,
        "mag": 0.0,
    },
]

params_shifts_list = [
    {
       "texture": "skx",
       # -- settings of periodic boundary conditions
       "system_sizes": [137],
       # -- general model parameters
       "t": -1,  # hopping
       "m": 5,  # exchange
       # -- phase shift diagram
       "q_shift": 24/137,  # fixed q value for phase shift diagram
       "n_shift": 300,  # number of phaseshifts for which the spectra are computed
       "min_shift": 0.25*pi,  # lower bound for range of shifts (included)
       "max_shift": 0.45*pi,  # upper bound for range of shifts (included)
       # -- net magnetization (up to normalisation factor)
       "mag": 0.7,
       # -- eta for regularisation of normalisation
       "eta": 0.0,
    },
]

params_shifts_kpm_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0/np.sqrt(3)),pi+np.arccos(0/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.7,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.7/np.sqrt(3)),pi+np.arccos(0.7/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # helical skyrmion lattice with regularisation
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0/np.sqrt(3)),pi+np.arccos(0/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.7,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.7/np.sqrt(3)),pi+np.arccos(0.7/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0),pi+np.arccos(0)],
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- eta for regularisation of normalisation
        "eta": 0.,
    },
    # sinusoidal skyrmion lattice with regularisation
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0),pi+np.arccos(0)],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 500,  # number of phaseshifts for which the spectra are computed
        "min_shift": 0,  # lower bound for range of shifts (included)
        "max_shift": pi,  # upper bound for range of shifts (included)
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
]

params_mags_list = [
]

params_mags_kpm_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 500,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(np.sqrt(3)*1/np.sqrt(3)),pi+np.arccos(np.sqrt(3)*1/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # helical skyrmion lattice with regularisation
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 500,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(np.sqrt(3)*1/np.sqrt(3)),pi+np.arccos(np.sqrt(3)*1/np.sqrt(3))],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "shift": 1/4*pi,
        "n_mag": 500,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # sinusoidal skyrmion lattice with regularisation
        {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 401,
        "q": 70/401,  # fixed q value for phase shift diagram
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,  # exchange
        # -- settings kpm
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        # -- positions in phase shift-magnetization phase diagram
        "shift": 1/4*pi,
        "n_mag": 500,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
]

params_chern_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": [57 for _ in range(6)],
        "q": [10/57 for _ in range(6)], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(6)],
        "mag": [0.0 for _ in range(6)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "fermi": [-9.94, -9.42, -1.21, 0.27, 0.92, 1.48], #theta: ~0.175
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": [57 for _ in range(8)],
        "q": [10/57, 10/57, 24/57, 10/57, 10/57, 10/57, 10/57, 24/57], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(8)],
        "mag": [0.0 for _ in range(8)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "fermi": [-9.46, -8.35, -4.55, -1.21, 0.90, 1.93, 2.88, 6.10], #theta: ~0.175 and 0.421
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    # vortex lattice
    {
        "texture": "sdw",
        # -- settings of periodic boundary conditions
        "system_sizes": [57 for _ in range(4)],
        "q": [10/57, 10/57, 10/57, 25/57], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(4)],
        "mag": [0.0 for _ in range(4)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "fermi": [-6.75, -5.50, -4.80, -2.03], #theta: ~0.175
        # -- delta for difference quotient
        "delta": 10**-5,
    },
]

params_chern_ids_list = [
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": [47 for _ in range(1,16)],
        "q": [(8,47) for _ in range(1,16)], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(1,16)],
        "mag": [0.0 for _ in range(1,16)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "ids": [1+(p/q)**2+0.5/q**2 for p,q in [(8,47) for _ in range(1,16)] ],
        # -- delta for difference quotient
        "delta": [10**-i for i in range(1,16)],
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": [59 for _ in range(1,16)],
        "q": [(10,59) for _ in range(1,16)], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(1,16)],
        "mag": [0.0 for _ in range(1,16)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "ids": [1+(p/q)**2+0.5/q**2 for p,q in [(10,59) for _ in range(1,16)] ],
        # -- delta for difference quotient
        "delta": [10**-i for i in range(1,16)],
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": [17,23,29,41,47,53,59,71,83,89,101,107,113],
        "q": [(3,17),(4,23),(5,29),(7,41),(8,47),(9,53),(10,59),(12,71),(14,83),(15,89),(17,101),(18,107),(19,113)], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(17)],
        "mag": [0.0 for _ in range(17)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "ids": [1+(p/q)**2+0.5/q**2 for p,q in [(3,17),(4,23),(5,29),(7,41),(8,47),(9,53),(10,59),(12,71),(14,83),(15,89),(17,101),(18,107),(19,113)]],
        # -- delta for difference quotient
        "delta": [10**-5 for _ in range(17)],
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": [17,23,29,41,47,53,59,71,83,89,101,107,113],
        "q": [(3,17),(4,23),(5,29),(7,41),(8,47),(9,53),(10,59),(12,71),(14,83),(15,89),(17,101),(18,107),(19,113)], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        # -- position in phase shift-magnetization phase diagram
        "shift": [0.0 for _ in range(17)],
        "mag": [0.0 for _ in range(17)],
        # -- set of labels from whose subsets of even cardinality the Chern numbers are calculated
        "tau1": True,
        "tau2": True,
        "u1": True,
        "u2": True,
        "u3": False,
        # -- fermi energy
        "ids": [1+(p/q)**2+0.5/q**2 for p,q in [(3,17),(4,23),(5,29),(7,41),(8,47),(9,53),(10,59),(12,71),(14,83),(15,89),(17,101),(18,107),(19,113)]],
        # -- delta for difference quotient
        "delta": [10**-8 for _ in range(17)],
    },
]

params_chern_fermis_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 22/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 53,
        "q": 20/53, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.7,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 53,
        "q": 20/53, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.7,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 24/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 53,
        "q": 20/53, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": pi,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": -0.5,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 53,
        "q": 20/53, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": -0.5,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-8,
    },
    # vortex lattice
    {
        "texture": "sdw",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -9.0,
        "max_fermi": 6.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "sdw",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 25/57, # fixed q value for Chern number calculation
        "flux": 0.0, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 5,   # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -9.0,
        "max_fermi": 6.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    # ferromagnet with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 0.0, # fixed q value for Chern number calculation
        "flux": 10/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    # skyrmion lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 10/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    # meron-antimeron lattice with finite flux
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 10/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": pi/2,
        "mag": 0.0,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "tau1u1": True, 
        "tau1u2": True,
        "tau2u1": True,
        "tau2u2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -11.0,
        "max_fermi": 8.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
]

params_chern_shifts_ids_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0/np.sqrt(3)),pi+np.arccos(0/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.7,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.7/np.sqrt(3)),pi+np.arccos(0.7/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0,
    },
    # helical skyrmion lattice with regularisation
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0/np.sqrt(3)),pi+np.arccos(0/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.7,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.7/np.sqrt(3)),pi+np.arccos(0.7/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0),pi+np.arccos(0)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # sinusoidal skyrmion lattice with regularisation
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0),pi+np.arccos(0)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 1000,
        "min_shift": 0.0,
        "max_shift": pi,
        "mag": 0.25,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.25),pi+np.arccos(0.25)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
]

params_chern_mags_ids_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(np.sqrt(3)*1/np.sqrt(3)),pi+np.arccos(np.sqrt(3)*1/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0,
    },
    # helical skyrmion lattice with regularisation
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": 0.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(np.sqrt(3)*1/np.sqrt(3)),pi+np.arccos(np.sqrt(3)*1/np.sqrt(3))],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0.1,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": -1.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": -1.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # sinusoidal skyrmion lattice with regularisation
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": -1.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "shift": 0.0,
        "n_mag": 1000,
        "min_mag": -1.0,
        "max_mag": 1.0,
        # -- origin in phase space
        "phi": [pi+np.arccos(0.5),pi+np.arccos(0.5)],
        # -- save spectrum for every point
        "save_spec": True,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
]

params_chern_diagram_list = [
    # helical skyrmion lattice
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 100,
        "min_shift": 0.0,
        "max_shift": pi,
        "n_mag": 100,
        "min_mag": 0.0,
        "max_mag": 2.0,
        # -- origin in phase space
        "phi": [0.0,0.0],
        # -- save spectrum for every point
        "save_spec": False,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 100,
        "min_shift": 0.0,
        "max_shift": pi,
        "n_mag": 100,
        "min_mag": 0.0,
        "max_mag": 2.0,
        # -- origin in phase space
        "phi": [0.0,0.0],
        # -- save spectrum for every point
        "save_spec": False,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+3*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-8,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    # sinusoidal skyrmion lattice
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 100,
        "min_shift": 0.0,
        "max_shift": pi,
        "n_mag": 100,
        "min_mag": 0.0,
        "max_mag": 2/np.sqrt(3),
        # -- origin in phase space
        "phi": [0.0,0.0],
        # -- save spectrum for every point
        "save_spec": False,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+1*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
    {
        "texture": "sin",
        # -- settings of periodic boundary conditions
        "system_sizes": 57,
        "q": 10/57, # fixed q value for Chern number calculation
        "flux": 0/57, # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        # -- positions in phase shift-magnetization phase diagram
        "n_shift": 100,
        "min_shift": 0.0,
        "max_shift": pi,
        "n_mag": 100,
        "min_mag": 0.0,
        "max_mag": 2/np.sqrt(3),
        # -- origin in phase space
        "phi": [0.0,0.0],
        # -- save spectrum for every point
        "save_spec": False,
        # -- set of Chern numbers to be calculated
        "tau1tau2": True,
        "u1u2": True,
        "tau1tau2u1u2": True,
        # -- orbital magnetism
        "morb": True,
        # -- IDS of gap
        "ids": 1+2*(10/57)**2+0.5/57**2,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
]

params_morb_scaling_list = [
    {
        "texture": "skx",
        # -- settings of periodic boundary conditions
        "system_sizes": 59,
        "q": [i/59 for i in range(0,29,1)], # fixed q value for Chern number calculation
        "flux": [i/59 for i in range(0,29,2)], # fixed flux value for Chern number calculation
        # -- general model parameters
        "t": -1,     # hopping
        "m": 5,      # exchange
        "kbT": 0, # finite temperature
        # -- position in phase shift-magnetization phase diagram
        "shift": 0.0,
        "mag": 0.0,
        # -- ids gap labels
        "gap_pos": "middle",
        "block": 1,
        "tau1tau2": 1,
        "tau1tau2u1u2": -1,
        # -- delta for difference quotient
        "delta": 10**-5,
        # -- eta for regularisation of normalisation
        "eta": 0.0,
    },
]
