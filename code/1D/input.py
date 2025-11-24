# ----------------------------------------------------------------------
# setup model
# ----------------------------------------------------------------------

params_theta_list = [
    {
        "texture": "limacon",
        # -- settings of periodic boundary conditions
        "system_sizes": [6421],
        # -- general model parameters
        "t": -1,  # hopping
        "m": 1,   # exchange
    },   
]

params_theta_kpm_list = [
    {
        "texture": "limacon",
        # -- settings of periodic boundary conditions
        "system_sizes": 160000,
        "n_q": 400,
        # -- settings kpm
        "scale": 4,
        "n_energies": 2048,
        "n_moments": 2048,
        "n_random_states": 10,
        "epsilon": 0.01,
        # -- general model parameters
        "t": -1,  # hopping
        "m": 1,  # exchange
    },   
]

params_chern_list = [
    {
        "texture": "limacon",
        # -- settings of periodic boundary conditions
        "system_sizes": [3249,3249,3249,3249,3249,3249],
        "q": [357/3249,357/3249,357/3249,835/3249,835/3249,835/3249], # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 1,   # exchange
        "kbT": 0, # finite temperature
        # -- fermi energy
        "fermi": [-3.2,0.0,3.2,-1.0,0.0,1.0],
        # -- delta for difference quotient
        "delta": 10**-5,
    },   
]

params_fermi_list = [
    {
        "texture": "limacon",
        # -- settings of periodic boundary conditions
        "system_sizes": 3249,
        "q": 357/3249, # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 1,   # exchange
        "kbT": 0, # finite temperature
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -4.0,
        "max_fermi": 4.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
    {
        "texture": "limacon",
        # -- settings of periodic boundary conditions
        "system_sizes": 3249,
        "q": 835/3249, # fixed q value for Chern number calculation
        # -- general model parameters
        "t": -1,  # hopping
        "m": 1,   # exchange
        "kbT": 0, # finite temperature
        # -- fermi energy
        "n_fermi": 2000,
        "min_fermi": -4.0,
        "max_fermi": 4.0,
        # -- delta for difference quotient
        "delta": 10**-5,
    },
]