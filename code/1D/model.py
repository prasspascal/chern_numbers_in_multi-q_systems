import numpy as np
from numba import jit
from constants import pi,sx,sy,sz


#-- hilbert space -----------------------------------------------------

@jit(nopython=True)
def construct_hilbert_space(n_x):
   """
      systematic labeling of quantum states
   """

   labels = np.zeros( (n_x*2, 2), dtype=np.integer )
   states = np.zeros( (n_x, 2), dtype=np.integer  )
   
   #construct hilbert space
   ii = 0
   for i in range(n_x):
      for s in range(2): #spin
         labels[ii]  = np.array( [i,s] )
         states[i,s] = ii
         ii += 1

   return labels, states


#-- magnetic texture --------------------------------------------------

@jit(nopython=True)
def magnetization(i, theta_1, phi_1, theta_2, phi_2, gamma=1.0, texture='limacon'):
   """
      Unraveling of the noncommutative torus onto the skyrmion lattice
   """

   n = np.zeros(3, dtype=np.float64)

   if texture=='limacon':
      # spin chain double-Q phase
      c1 = np.cos( 2*pi * ( i* theta_1 ) + phi_1 )
      s1 = np.sin( 2*pi * ( i* theta_1 ) + phi_1 )

      c2 = np.cos( 2*pi * ( i* theta_2 ) + phi_2  )
      s2 = np.sin( 2*pi * ( i* theta_2 ) + phi_2  )

      n = np.array([0,-s1,c1], dtype=np.float64) + gamma * np.array([0,-s2,c2], dtype=np.float64)

   elif texture=='fm':
      n = np.array([0,0,1], dtype=np.float64)

   return n


#-- hamiltonian setup -------------------------------------------------

@jit(nopython=True)
def nearest_neighbors(i, n_x, pbc):

   nn = []

   if pbc:
      nn.append( (i+1) % n_x )
   elif (i+1)<n_x: 
      nn.append( (i+1) )

   return np.array(nn)

@jit(nopython=True)
def set_hamiltonian(hamiltonian, states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc):
   """
      nearest neighbor hopping plus exchange
   """

   for i in range(n_x):
      #hopping 
      nn = nearest_neighbors(i, n_x, pbc)
      for n in nn:
         for s in range(2):
            ii = states[i,s]
            jj = states[n,s]

            hamiltonian[ii, jj] = t
            hamiltonian[jj, ii] = t

      #onsite term
      nvec   = magnetization(i, theta_1, phi_1, theta_2, phi_2, gamma, texture)
      onsite = m * (nvec[0]*sx + nvec[1]*sy + nvec[2]*sz)
      for a in range(2):
         for b in range(2):
            ii = states[i,a]
            jj = states[i,b]
            hamiltonian[ii, jj] = onsite[a,b]

def set_hamiltonian_kpm(hamiltonian, states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc):
   """
      nearest neighbor hopping plus exchange
   """

   for i in range(n_x):
      #hopping 
      nn = nearest_neighbors(i, n_x, pbc)
      for n in nn:
         for s in range(2):
            ii = states[i,s]
            jj = states[n,s]

            hamiltonian[ii, jj] = t
            hamiltonian[jj, ii] = t

      #onsite term
      nvec   = magnetization(i, theta_1, phi_1, theta_2, phi_2, gamma, texture)
      onsite = m * (nvec[0]*sx + nvec[1]*sy + nvec[2]*sz)
      for a in range(2):
         for b in range(2):
            ii = states[i,a]
            jj = states[i,b]
            hamiltonian[ii, jj] = onsite[a,b]


#-- linear algebra ----------------------------------------------------

def spectrum(states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc):
   """
      calculate the spectrum by exact diagonalization
   """
   
   dim = 2*n_x
   H = np.zeros((dim,dim), dtype=np.complex128)
   set_hamiltonian(H, states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc)

   return  np.linalg.eigvalsh(H, UPLO='U').real

def spectrumV(states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc):
    """
    calculate the spectrum by exact diagonalization, return eigenvalues and eigenvectors
    """

    dim = 2*n_x
    H = np.zeros((dim, dim), dtype=np.complex128)
    set_hamiltonian(H, states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc)

    return np.linalg.eigh(H, UPLO="U")

@jit(nopython=True)
def spec_projection(eigvals, eigvecs, fermi,kbT=0):
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
def fermi_distribution(E,mu,kbT):
    return 1.0 / ( np.exp( (E-mu) / kbT) + 1)

# -- derivations ----------------------------------------------------

@jit(nopython=True)
def translation_derivative(P, states, n_x):
   """
   Calculate the position derivation of the spectral projection

   P: spectral projection
   states: returns the index of the state with a specific label
   n_x: size of the unit cell
   """

   dP=np.zeros_like(P)

   for i in range(n_x):
      for s in range(2):
         for j in range(n_x):
            for t in range(2):

               row = states[i, s]
               col = states[j, t]

               dP[row, col] = P[row, col] * ((i-j)-int((i-j)/(0.5*n_x))*n_x)

   return -1j*dP

def phason_derivative(P_phi, dphi, states, n_x, theta_1, phi_1, theta_2, phi_2, gamma, t, m, texture, pbc, fermi,kbT=0):
   """
   Return the difference quotient of the spectral projection at phi and phi+dphi

   P       : spectral projection at phi
   dphi    : step width for directional derivative
   others  : parameters of spectrumV of P_phi
   """

   eigvals, eigvecs = spectrumV(states, n_x, theta_1, phi_1+dphi, theta_2, phi_2+theta_2/theta_1*dphi, gamma, t, m, texture, pbc)

   P_phidphi = spec_projection(eigvals, eigvecs, fermi, kbT)

   return ( P_phidphi - P_phi )/dphi

# -- trace ----------------------------------------------------

@jit(nopython=True)
def trace(V,A):
    """
    Compute trace of A
    """

    return 1/V * np.trace(A)

# -- chern character ------------------------------------------

@jit(nopython=True)
def chern_character(P,dP, n_x):
   """
   Compute the chern character for all subsets of even cardinality of the set of labels True in J
   """

   # dictionary for chern numbers
   chern_numbers = {}

   # always calculate the 0th chern number
   chern_numbers[""] = trace(n_x,P)

   # compute the chern number by computing twice the imaginary part of the trace of the product of derivations
   chern_numbers["tau1u1"] = (2*pi*1j) * 2*1j*trace(n_x, P.dot(dP[0].dot(dP[1]))).imag

   return chern_numbers