import numpy as np
from constants import pi, rot120, rot240, sx, sy, sz
from lattice import G, R, Twist, Twist3D
from numba import jit

# -- hilbert space -----------------------------------------------------

@jit(nopython=True)
def construct_hilbert_space(n_x, n_y):
    """
    systematic labeling of quantum states
    """

    labels = np.zeros((n_x * n_y * 2, 3), dtype=np.integer)
    states = np.zeros((n_x, n_y, 2), dtype=np.integer)

    # construct hilbert space
    ii = 0
    for i in range(n_x):
        for j in range(n_y):
            for s in range(2):  # spin
                labels[ii] = np.array([i, j, s])
                states[i, j, s] = ii
                ii += 1

    return labels, states


# -- magnetic texture --------------------------------------------------

@jit(nopython=True)
def magnetization(i, j, theta, alpha=0.0, phi=np.array([0.0,0.0]), shift=0.0, mag=0.0, texture="skx", eta=0.0):
    """
    Unraveling of the noncommutative torus onto the skyrmion lattice

    theta: wavenumber
    alpha: twist angle
    phi  : origin in phason (used for phason derivatives)
    shift: phase shift (applied to phase 1), relative phase degree of freedom moving the origin
    mag: parameter to tune net magnetization (proportional up to normalization)
    """

    # real-space position
    x = R(i, j)

    # magnetization vector
    n = np.zeros(3, dtype=np.float64)

    # rotation of magnetic lattice w.r.t atomic lattice
    twist = Twist(alpha)
    twist3D = Twist3D(alpha)

    # phase factors
    phase_1 = theta * np.dot(x, np.dot(twist, G( 1,  0))) + phi[0] + shift
    phase_2 = theta * np.dot(x, np.dot(twist, G( 0,  1))) + phi[1]
    phase_3 = theta * np.dot(x, np.dot(twist, G(-1, -1))) - phi[0] - phi[1]

    if texture == "skx":
        # 'skyrmion' helical triple-Q phase

        n1 = np.dot(rot240, np.array([np.sin(phase_1), 0, np.cos(phase_1)]))
        n2 = np.array([np.sin(phase_2), 0, np.cos(phase_2)])
        n3 = np.dot(rot120, np.array([np.sin(phase_3), 0, np.cos(phase_3)]))

        n = np.dot(twist3D, n1 + n2 + n3) + np.array([0, 0, np.sqrt(3) * mag])
        n = n * 1 / (np.linalg.norm(n) + eta)

    elif texture == "sin":
        # sinusoidal triple-Q phase

        # Gamma = 0 and theta = arccos( 1/sqrt(3) )
        n1 = np.dot(rot240, np.array([-np.sqrt(2/3)*np.cos(phase_1), 0, np.sqrt(1/3)*np.cos(phase_1)]))
        n2 = np.array([-np.sqrt(2/3)*np.cos(phase_2), 0, np.sqrt(1/3)*np.cos(phase_2)])
        n3 = np.dot(rot120, np.array([-np.sqrt(2/3)*np.cos(phase_3), 0, np.sqrt(1/3)*np.cos(phase_3)]))

        n = np.dot(twist3D, n1 + n2 + n3) + np.array([0, 0, np.sqrt(3) * mag])
        n = n / (np.linalg.norm(n) + eta)

    elif texture == "sdw":
        # spin-density wave triple-Q phase

        n1 = np.dot(rot240, np.array([0, np.sin(phase_1), 0]))
        n2 = np.array([0, np.sin(phase_2), 0])
        n3 = np.dot(rot120, np.array([0, np.sin(phase_3), 0]))

        n = np.dot(twist3D, n1 + n2 + n3) / 3.0

    elif texture == "fm":
        # ferromagnet

        n = np.array([0, 0, 1], dtype=np.float64)

    return n


# -- hamiltonian setup -------------------------------------------------
@jit(nopython=True)
def nearest_neighbors(i, j, n_x, n_y, boundary_conditions):
    """
    Returns a list of nearest neighbors to to the site i,j.
    The last number in each entry specify which translation was performed
    """

    if boundary_conditions == "periodic":

        nn = [[(i + 1) % n_x, j], [i, (j + 1) % n_y], [(i + 1) % n_x, (j - 1) % n_y]]

    elif boundary_conditions == "open":

        nn = []

        # check if neighbors exist
        ni = i + 1
        if ni < n_x:
            nn.append([ni, j])

        nj = j + 1
        if nj < n_y:
            nn.append([i, nj])

        nj2 = j - 1
        if 0 <= nj2 < n_y and ni < n_x:
            nn.append([ni, nj2])

    else:
        raise ValueError(
            "nearest_neighbors: unkown boundary condition( choose either 'periodic' or 'open')"
        )

    return nn


@jit(nopython=True)
def set_hamiltonian(
    hamiltonian,
    states,
    n_x,
    n_y,
    theta,
    alpha,
    phi,
    shift,
    mag,
    t,
    m,
    texture,
    boundary_conditions,
    regularisation = 0.0,
    flux = 0.0
):
    """
    nearest neighbor hopping on the triangular lattice plus exchange
    """

    for i in range(n_x):
        for j in range(n_y):

            # hopping
            nn  = nearest_neighbors(i, j, n_x, n_y, boundary_conditions)
            hopp = [ t* np.exp(1j* pi* flux * -j), t* np.exp(1j* pi* flux * +i), t* np.exp(1j* pi* flux * -(i+j))]
            for s in range(2):
                ii = states[i, j, s]
                for site,h in zip(nn,hopp):
                    jj = states[site[0], site[1], s]

                    hamiltonian[ii, jj] = h
                    hamiltonian[jj, ii] = np.conjugate(h)

            # onsite term
            nvec = magnetization(i, j, theta, alpha, phi, shift, mag, texture, regularisation)
            onsite = m * (nvec[0] * sx + nvec[1] * sy + nvec[2] * sz)
            for a in range(2):
                for b in range(2):
                    ii = states[i, j, a]
                    jj = states[i, j, b]
                    hamiltonian[ii, jj] = onsite[a, b]
                    hamiltonian[jj, ii] = onsite[b, a]

def set_hamiltonian_kpm(
    hamiltonian,
    states,
    n_x,
    n_y,
    theta,
    alpha,
    phi,
    shift,
    mag,
    t,
    m,
    texture,
    boundary_conditions,
    regularisation = 0.0,
    flux = 0.0
):
    """
    nearest neighbor hopping on the triangular lattice plus exchange
    """

    for i in range(n_x):
        for j in range(n_y):

            # hopping
            nn  = nearest_neighbors(i, j, n_x, n_y, boundary_conditions)
            hopp = [ t* np.exp(1j* pi* flux * -j), t* np.exp(1j* pi* flux * +i), t* np.exp(1j* pi* flux * -(i+j))]
            for s in range(2):
                ii = states[i, j, s]
                for site,h in zip(nn,hopp):
                    jj = states[site[0], site[1], s]

                    hamiltonian[ii, jj] = h
                    hamiltonian[jj, ii] = np.conjugate(h)

            # onsite term
            nvec = magnetization(i, j, theta, alpha, phi, shift, mag, texture, regularisation)
            onsite = m * (nvec[0] * sx + nvec[1] * sy + nvec[2] * sz)
            for a in range(2):
                for b in range(2):
                    ii = states[i, j, a]
                    jj = states[i, j, b]
                    hamiltonian[ii, jj] = onsite[a, b]
                    hamiltonian[jj, ii] = onsite[b, a]


# -- linear algebra ----------------------------------------------------

def spectrum(
    states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, regularisation = 0.0, flux = 0.0
):
    """
    calculate the spectrum by exact diagonalization, return eigenvalues
    """

    dim = 2 * n_x * n_y
    H = np.zeros((dim, dim), dtype=np.complex128)
    set_hamiltonian(
        H, states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, regularisation, flux
    )

    return np.linalg.eigvalsh(H, UPLO="U").real

def spectrumV(
    states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, regularisation = 0.0, flux = 0.0
):
    """
    calculate the spectrum by exact diagonalization, return hamiltonian, eigenvalues, and eigenvectors
    """

    dim = 2 * n_x * n_y
    H = np.zeros((dim, dim), dtype=np.complex128)
    set_hamiltonian(
        H, states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, regularisation, flux
    )

    eigvals, eigvecs = np.linalg.eigh(H, UPLO="U")

    return H, eigvals, eigvecs

@jit(nopython=True)
def fermi_projection(
    eigvals, eigvecs, fermi, kbT=0
):
    """
    Create a projection D_P onto the occupied subspace in eigenvector basis and switch to standard basis
    P = V D_P V^H 
    with V   = (v1, ..., vN) matrix with normalized eigenvectors in columns
         D_P = Diag(1, ..., 1, 0, ..., 0) diagonal projection of rank equal to number of eigenvalues below fermi
         ^H  hermitian conjugate i.e. transpose and complex conjugate
         P   projection onto eigenspaces with eigenvalues below fermi
    """

    # array with 1 for each occupied state (energy less or equal fermi energy) and 0 for each unoccupied state (energy greater fermi energy)
    if kbT==0:
        occ = np.array([1 if eigval<= fermi else 0 for eigval in eigvals],dtype=np.complex128)
    else:
        occ = np.array([ fermi_distribution(eigval,fermi,kbT) for eigval in eigvals],dtype=np.complex128)

    # projection onto occupied subspace of states (eigenbasis)
    eig_projection = np.diag(occ)

    # return projection onto occupied subspace of states (statebasis)
    return eigvecs.dot(eig_projection.dot(eigvecs.transpose().conj()))

@jit(nopython=True)
def spec_projection(
    n_x, n_y, eigvecs, ids
):
    """
    Create a projection D_P onto the occupied subspace in eigenvector basis and switch to standard basis
    P = V D_P V^H
    with V   = (v1, ..., vN) matrix with normalized eigenvectors in columns
         D_P = Diag(1, ..., 1, 0, ..., 0) diagonal projection of rank equal to number of occupied
         ^H  hermitian conjugate i.e. transpose and complex conjugate
         P   projection onto occupied states in eigenbasis
    """

    # number of occupied states
    n_occ = int(ids*n_x*n_y)
    n_unocc = 2*n_x*n_y - n_occ

    # array with 1 for each occupied state and 0 for each unoccupied state
    occ = np.concatenate( (np.ones(n_occ,dtype=np.complex128),np.zeros(n_unocc,dtype=np.complex128)) )

    # projection onto occupied subspace of states (eigenbasis)
    eig_projection = np.diag(occ)

    # return projection onto occupied subspace of states (statebasis)
    return eigvecs.dot(eig_projection.dot(eigvecs.transpose().conj()))

@jit(nopython=True)
def fermi_distribution(E,mu,kbT):
    return 1.0 / ( np.exp( (E-mu) / kbT) + 1)

# -- derivations ----------------------------------------------------

@jit(nopython=True)
def translation_derivative(
    P, index, states, n_x, n_y
):
    """
    Calculate the position derivation of the spectral projection

    P: spectral projection
    index: index of partial derivative
    states: returns the index of the state with a specific label
    n_x, n_y: size of the unit cell
    """

    dP=np.zeros_like(P)
    #L_x = (n_x-1)/2
    #L_y = (n_y-1)/2

    for i in range(n_x):
        for j in range(n_y):
                for s in range(2):
                        for k in range(n_x):
                            for m in range(n_y):
                                for t in range(2):

                                    row = states[i, j, s]
                                    col = states[k, m, t]

                                    if index==0:
                                        #dP[row, col] = P[row, col] * ((i-k+L_x)%n_x-L_x)
                                        dP[row, col] = P[row, col] * ((i-k)-int((i-k)/(0.5*n_x))*n_x)
                                    elif index==1:
                                        #dP[row, col] = P[row, col] * ((j-l+L_y)%n_y-L_y)
                                        dP[row, col] = P[row, col] * ((j-m)-int((j-m)/(0.5*n_y))*n_y)
    return -1j*dP

def phason_derivative_DQ_fermi(
    P_phi, dphi, states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, fermi, kbT=0, regularisation=0.0, flux = 0.0
):
    """
    Return the difference quotient of the fermi projection at phi and phi+dphi
    IMPORTANT (at zero temperature):
    The fermi energy for a fixed filling is not the same between the spectral projection at phi and phi+dphi and needs to be reevaluated accordingly to get projectors of the same rank

    P       : fermi projection at phi
    dphi    : array of step width for directional derivative
    others  : parameters of spectrumV of P_phi
    """

    # diagonalise displaced hamiltonian
    _, eigvals, eigvecs = spectrumV(states, n_x, n_y, theta, alpha, phi+dphi, shift, mag, t, m, texture, boundary_conditions, regularisation, flux)

    if kbT==0:
        # fermi energy of occupied states with filling as in P_phi (this ensures the same rank between P_phi and P_phidphi)
        filling = np.trace(P_phi)
        fermi = eigvals[int(filling)-1]

    # compute displaced projection
    P_phidphi = fermi_projection(eigvals, eigvecs, fermi, kbT)

    return ( P_phidphi - P_phi )/np.linalg.norm(np.array(dphi))

def phason_derivative_DQ_IDS(
    P_phi, dphi, states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions, ids, regularisation=0.0, flux = 0.0
):
    """
    Return the difference quotient of the spectral projection at phi and phi+dphi

    P       : spectral projection at phi
    dphi    : array of step width for directional derivative
    ids     : integrated density of states the spectral projection projects to
    others  : parameters of spectrumV of P_phi
    """

    # diagonalise displaced hamiltonian
    _, _, eigvecs = spectrumV(states, n_x, n_y, theta, alpha, phi+dphi, shift, mag, t, m, texture, boundary_conditions, regularisation, flux)

    P_phidphi = spec_projection(n_x, n_y, eigvecs, ids)

    return ( P_phidphi - P_phi )/np.linalg.norm(np.array(dphi))

# -- trace ----------------------------------------------------

@jit(nopython=True)
def trace(
        V,A
):
    """
    Compute trace of A
    """

    return 1/V * np.trace(A)
    

# -- chern character ------------------------------------------

@jit(nopython=True)
def chern_character(
        J,P,dP_J, n_x, n_y
):
    """
    Compute the chern character for all subsets of even cardinality of the set of labels True in J
    """

    # names of the labels
    labels = ["tau1","tau2","u1","u2","u3"]

    # dictionary for chern numbers
    chern_numbers = {}

    # always calculate the 0th chern number
    chern_numbers[""] = trace(n_x*n_y,P)

    # if there are at least 2 labels True in J, calculate all corresponding 1st Chern numbers
    if np.sum(J,dtype=np.int16)>=2:
        # for all combinations of 2 elements from the list of labels
        for a,b in comb([0,1,2,3,4],2):
            # check if both labels are True
            if J[a] and J[b]:
                # compute the chern number by computing the trace of all permutations of derivations corresponding to two True labels
                chern_numbers[labels[a]+labels[b]] = (2*pi*1j) * cyclic_cocycle(n_x*n_y,P,dP_J,[a,b])

    # if there are at least 4 labels True in J, calculate all corresponding 2nd Chern numbers
    if np.sum(J,dtype=np.int16)>=4:
        # for all combinations of 4 elements from the list of labels
        for a,b,c,d in comb([0,1,2,3,4],4):
            # check if both labels are True
            if J[a] and J[b] and J[c] and J[d]:
                # compute the chern number by computing the trace of all permutations of derivations corresponding to four True labels
                chern_numbers[labels[a]+labels[b]+labels[c]+labels[d]] = (2*pi*1j)**2 / 2 * cyclic_cocycle(n_x*n_y,P,dP_J,[a,b,c,d])

    return chern_numbers

@jit(nopython=True)
def cyclic_cocycle(V,P,dP,J):
    """
    Compute the cyclic cocycle as the trace of the product of the spectral projection and the sum over the product of the permutations of its derivatives (listed in label set J)
    """

    # variable for the product of the spectral projection and the sum over the product of the permutations of its derivatives (listed in label set J)
    A = np.zeros_like(P)

    # generate array of al permutations of labels in J
    permutations = [p for p in perm(J)]

    # remove tuples which are the reverse of another tuple in the list
    # the reverse permutation corresponds to the complex conjugate
    # if the sign of the reverse permutation is equal we take twice the real part as result
    # if the sign of the reverse permutation is switched we take twice the imaginary part as result
    permutations = remove_reverse_arrays(permutations)

    # iterate over permutations of labels in J (up to reverse order)
    for p in permutations:

        # variable for product of the derivatives (in order of one permutation p of labels in J)
        temp = np.zeros_like(dP[p[0]])
        temp += dP[p[0]]

        # compute product of the derivatives
        for i in range(1,len(p)):
            temp = temp.dot(dP[p[i]])

        temp *= perm_sign(p)

        # add the product to the sum
        A += temp

    # product with the spectral projection
    A = P.dot(A)

    # return twice the real or imaginary part of the trace (depending in the sign of the reverse permutation)
    if perm_sign(np.array(J)[::-1])==1:
        return 2*trace(V,A).real
    elif perm_sign(np.array(J)[::-1])==-1:
        return 2*1j*trace(V,A).imag
  
# -- orbital magnetisation ------------------------------------

@jit(nopython=True)
def orbital_magnetisation(
        H, P, dP, fermi, n_x, n_y
):
    """
    Compute the z-component of the orbital magnetisation contributions
    (up to prefactor e/(\hbar c) ) from modern theory

    return tuple: Morb, MorbLC, MorbIC, Ch_tau1tau2

    Morb: orbital magnetisation = MorbLC + MorbIC - 2 * fermi * Ch_tau1tau2
    MorbLC: local circulation contribution
    MorbIC: itinerant current contribution
    Ch_tau1tau2: momentum-space Chern number
    """

    # compute products
    dP1dP2 = dP[0].dot(dP[1])
    PdP1dP2 = P.dot(dP1dP2)

    # compute Chern number (2 pi i) (2i) Im Tr(P dP1 dP2)
    Ch_tau1tau2 = (-4*pi) * trace(n_x*n_y, PdP1dP2).imag

    # compute itinerant current contribution Im Tr(H P dP1 dP2)
    MorbIC = trace(n_x*n_y, H.dot(PdP1dP2)).imag

    # compute local circulation contribution Im Tr(P dP1 H dP2)
    # Im Tr(H P dP1 dP2 - H dP1 dP2) = - Im Tr(H (1-P) dP1 dP2) = - Im Tr(H dP1 P dP2) = - Im Tr(P dP2 H dP1) = Im Tr(P dP1 H dP2)
    MorbLC = MorbIC-trace(n_x*n_y, H.dot(dP1dP2)).imag 

    # compute orbital magnetisation (z-component)
    Morb = MorbLC + MorbIC - 2 * fermi * Ch_tau1tau2/(-4*pi)

    return (Morb,MorbLC,MorbIC,Ch_tau1tau2)

# -- tools ----------------------------------------------------

@jit(nopython=True)
def comb(pool, r):
    '''
    from https://docs.python.org/3/library/itertools.html#itertools.combinations
    modified for compatibility with numba
    '''

    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield [pool[i] for i in indices]
    while True:
        for i in range(r-1,-1,-1):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield [pool[i] for i in indices]

@jit(nopython=True)
def perm(pool, r=None):
    '''
    https://docs.python.org/3/library/itertools.html#itertools.permutations
    modified for compatibility with numba
    '''

    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
    yield [pool[i] for i in indices[:r]]
    while n:
        for i in range(r-1,-1,-1):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield [pool[i] for i in indices[:r]]
                break
        else:
            return

@jit(nopython=True)
def remove_reverse_arrays(array_list):
    """
    Remove arrays which are the reverse of another array already in the list
    """
    array_list = np.array(array_list)
    result = []
    for i in range(len(array_list)):
        is_reverse = False
        for j in range(0, i):
            if np.all(array_list[i] == array_list[j][::-1]):
                is_reverse = True
                break
        if not is_reverse:
            result.append(array_list[i])
    return result

@jit(nopython=True)
def perm_sign(perm):
    """
    Written by OpenAI:
    The idea is that if the permutation is already sorted, then no swaps are needed and the sign is 1.
    If a swap is needed to sort the permutation, the sign will be -1. If two or more swaps are needed, the sign will be 1 again, and so on.

    This implementation has time complexity of O(n log n) in average and worst case, because it uses bubble sort internally.
    """
    sign = 1
    n = len(perm)
    for i in range(n-1):
        for j in range(i+1, n):
            if perm[i] > perm[j]:
                perm[i], perm[j] = perm[j], perm[i]
                sign *= -1
    return sign

@jit(nopython=True)
def mysum(xs):
    """
    Sum function for numpy arrays compatible with numba
    """
    tmp = np.zeros_like(xs[0])
    for x in xs:
        tmp += x

    return tmp

# -- tests ----------------------------------------------------

if __name__ == "__main__":

    n_x = 10
    n_y = 10

    labels, states = construct_hilbert_space(n_x, n_y)

    # print("labels: ", labels)

    theta = 1
    alpha = 1
    phi = [0,0,0]
    shift = 0
    mag = 0
    t = 1
    m = 1
    texture = "skx"
    boundary_conditions = "open"

    spec = spectrum(
        states, n_x, n_y, theta, alpha, phi, shift, mag, t, m, texture, boundary_conditions
    )

    print(np.amin(spec), np.amax(spec))
