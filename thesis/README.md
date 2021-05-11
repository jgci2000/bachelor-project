# Thesis Outline (30 pages)
-----------------------------------------

## Chap. 1 - Introduction (/ pages)
* State of Art 
* Introduce the problem
* Thesis outline
* Why is this important?
-----------------------------------------

## Chapter 2 - Spin Models (/ pages)

### Section 2.1 - Ising Model
* Introduce the Ising Model

#### 2.1.1 - Joint Density of States
 * What is the DOS and JDOS?

#### 2.1.2 - Thermodynamics
 * Thermodynamic relations computed from the JDOS

#### 2.1.3 - Relevance
 * Where is the Model applied?
   * Physics/Material Sciences -> Ferromagnetism
   * Ising Social
   * Finances

### Section 2.2 - Spin S Ising Model
* Introduce the Model
* Differences from Ising Model
-----------------------------------------

## Chapter 3 - Monte-Carlo Methods for Spin Models (/ pages)

### Section 3.1 - Metropolis

### Section 3.2 - Random Path Sampling

### Section 3.3 - Wang-Landau

-----------------------------------------
## Chapter 4 - Introduction to the Flat Scan Sampling (/ pages)

### Section 4.1 - Background (1/2 pages)
* How the method was created
* Its relation to RPS and WL
#### 4.1.1 - Algorithm
* How it actually works

### Section 4.2 - Implementation for the Ising model (/ pages)

#### 4.2.1 - Single Core
* Single Core C++
* Efforts to make the simulation more accurate -> skip

#### 4.2.2 - MPI
* One manager and N-1 workers
* One manager and N workers
* Efforts to make the simulation more accurate -> WL_shuffle

#### 4.2.3 - Performance and Scalling
* Comparison of Wall times for different lattices and sizes (for all implementations, MATLAB inc.)
* Multi-Core scalling of the method
-----------------------------------------

## Chapter 5 - Validation and Convergence of the FSS (/ pages)

### Section 5.1 - 
* Convergence
* Validation
* Plots of statistical analysis and REP analysis
* How does it compare to WL (comparison with WL)








## Chapter 5 - Finite Size Scalling and Thermodynamics (4/5 pages)


## Chapter 6 - Extension to the Ising SpinS Model (5/6 pages)
* Explain the different approach to the problem
* Implementation details
* Results (compare with JA Metropolis results)

## Chapter 7 - Conclusion and Future Work (1 page)












