# Non-commutative K-Theory of Lattice-Aperiodic Multi-q Magnetic Systems

### Code Base for the PhD Thesis by Pascal Jakob Alexander Praß

Date of defence: 29th of August 2025

## Abstract
Non-collinear magnetic structures provide a promising platform for energy-efficient carriers of information.
Periodic non-collinear spin arrangements generated as the interference pattern of spin waves are an extensive subclass known as multi-$q$ magnets.
The study of their electronic properties is of central importance for the implementation in future technology.
As the length scale of a multi-$q$ spin texture approaches the lattice constant of its host material, gapped topological states may form in the associated electronic system, similar to the formation of Landau levels in the presence of magnetic fields.
Given that the textures at these length scales are discrete and lattice incommensurate, we challenge the prevailing notion that their description relies on emergent magnetic fields.
Instead, we adopt a $C^\ast$-algebraic viewpoint that harmonises with these properties and facilitates the computation of invariants associated with the topology of the electronic states.
We implement a computational programme of non-commutative $K$-theory to compute all Chern numbers associated with real, momentum and mixed space.
As a central application, we tune texture parameters to create discontinuous jumps in the real space winding number of skyrmionic textures and observe the relation to the Chern numbers.
We find a peculiar discrepancy between the behaviour of the momentum space Chern number and the winding number, and identify a single energy level that contributes to the real space Chern number.
The non-commutative framework also provides access to the computation of the orbital magnetisation in accordance with modern theory even in the presence of finite magnetic flux.
This allows us to compare the orbital magnetisation with the texture and electronic state topology and investigate the scaling behaviour with respect to the texture length scale and magnetic flux.

## About this repository

The repository is divided into two parts. 

1. `/code` contains the Python source code which sets up the spin Hamiltonian, constructs the Fermi projections, and evaluates the Chern numbers and orbital magnetisation.
2. `/gfx` contains both the data analysis and the data sets that went into the creation of the figures in the manuscript.