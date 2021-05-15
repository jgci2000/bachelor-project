# Thesis Outline (30 pages)
-----------------------------------------

## Chap. 1 - Introduction (3/4 pages)
* State of Art 
* Introduce the problem
* Thesis outline
* Why is this important?
-----------------------------------------

## Chapter 2 - Spin Models (3/5 pages)

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

## Chapter 3 - Monte-Carlo Methods for Spin Models (4/5 pages)
* Describe the methods (...)

### Section 3.1 - Metropolis

### Section 3.2 - Random Path Sampling

### Section 3.3 - Wang-Landau

-----------------------------------------
## Chapter 4 - Introduction to the Flat Scan Sampling (5/6 pages)

### Section 4.1 - Background
* How the method was created
* Its relation to RPS and WL
#### 4.1.1 - Algorithm
* How it actually works

### Section 4.2 - Implementation for the Ising model
#### 4.2.1 - Single Core
* Single Core C++
* Efforts to make the simulation more accurate -> skip
#### 4.2.2 - MPI
* One manager and N-1 workers
* One manager and N workers
* Efforts to make the simulation more accurate -> WL_shuffle
#### 4.2.3 - Performance and Scaling
* Comparison of Wall times for different lattices and sizes (for all implementations, MATLAB inc.)
* Multi-Core scalling of the method
* Talk about Amdhal's Law and how it fits to the results
-----------------------------------------

## Chapter 5 - Validation and Convergence of the FSS (4/5 pages)
* Convergence
* Validation
* Plots of statistical analysis and REP analysis

### Section 5.1 - Convergence

### Section 5.2 - Validation
------------------------------------------

## Chapter 6 - Wang Landau Comparison (3/2 pages)
* Why compare with wang landau?
* good and bad for both methods
* How does it compare to WL (comparison with WL)
------------------------------------------

## Chapter 7 - Finite Size Scalling and Thermodynamics (/ pages)
------------------------------------------

## Chapter 8 - Extension to the Ising SpinS Model (/ pages)
* Explain the different approach to the problem
* Implementation details
* Results (compare with JA Metropolis results)
------------------------------------------

## Chapter 9 - Conclusion and Future Work (1 page)












