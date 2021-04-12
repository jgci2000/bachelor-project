# Project

This is my 3rd year Research Project in Monte-Carlo methods in Statistical Mechanics. It means that I use Monte-Carlo methods to solve models like the Ising Model, Heisenberg Model, ..., and get the Joint Density of States.

## Methods

* Wang Landau
* Metropolis
* Random Path Sampling
* Flat Scan Sampling

# Objectives
## Week 1:
* Compile code at the HCP
* Debug the FSS C++ code for large systems (check the histogram)
* Install MATLAB of the ryzen JA
* Run FSS MATLAB code for 3D SC L8 system

## Week 2: 
* Debug the FSS Ising-1/2 code
* Submit a job in the cluster

## Week 5:
* Optimize new FSS Ising-1/2 implementation
* Implement MPI version of my new FSS Ising-1/2 implementation
* Run simulations to validate my implementation
* Make a jupyter notebook for last point

## Week 6:
* Refurbish WL C++ implementation
* Run WL and FSS L16_SS to compare
* Continue developing REP and stat analysis scripts
* Benchmark and validade computations in the cluster
* Explore shuffle and skip effects trough 10 simulations of each (L16_SS R1E4)

## Week 7:
* makefiles and generate folder to write results to
* Script for shuffle_skip analysis
* Study the effects of the shuffle and skip


# TODO
## Important

* What should I change in a parallel FSS simulation so the results are the same as a single core simulation with the same parameters
  * Do a WL random walk before beggining the method (shuffle)
  * Change the skip parameter
  * Compare the results with 10-100 single core simulations
* Cluster
  * Benchmark node vs multi-node
  * Validate computations in the cluster, single core and parallel
* REP analysis
* skip analysis
* Is it worth to do a WL step before the FSS simulation to know all of the points in the phase space?
* WL vs FSS script

## Secondary

* Compare FSS with WL
  * Update the C++ WL implementation
  * Run simulations for L16_SS with various REP and f_final values
  * Compare var, var(1/sqrt(...)), runtime, ...
* Debug Ising Spin S code
  * functions working: scan_at_q1, random_config_at_q, WK_rw_criteria; where the problem could be -> scan or rw_step_at_q
* Thermodynamics and Finite Size effect scripts to Jupyter Notebooks
* Study the thermodynamics from different lattices
* Learn and study the critical exponents from lattices

