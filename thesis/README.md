# Report Outline (30 pages)
-----------------------------------------

## Chapter 1 - Introduction (3/4 pages)

* Computations in physics 
* Metropolis Method and why it is now perfect
* Efforts in the 90s to try and create new and faster methods
* Wang-Landau sampling and how it solved some problems
* Its applications and more importantly applications to ferromagnets (briefly explaing ferromagnets)
* Gaps that the WL did not fulfil
* Introduce the FSS 

* Report outline

## Plots for Chapter 1
* Ferromagnet plot from JA ppt (magnetization curve)
* Maybe a plot of Metropolis results compared with experimental results
* Plot of the DOS computations from the original paper

-----------------------------------------

## Chapter 2 - Ising Model (3/5 pages)
* Introduce the Ising Model
* History of the model
* Phase transition at Tc
* Ferromagnetism and Paramagnetism states of the Ising Model

### Section 2.1 - Joint Density of States
* What is the DOS and JDOS?
### Section 2.2 - Thermodynamics
* Thermodynamic relations computed from the JDOS
### Section 2.3 - Relevance
* Where is the Model applied?
  * Physics/Material Sciences -> Ferromagnetism
  * Ising Social
  * Finances

## Plots for Chapter 2
* Ferromagnet plot form JA ppt (oriented spins)
* JDOS table or the L2 system
* JDOS_exact plotted for L4

-----------------------------------------

## Chapter 3 - Monte-Carlo Methods for Spin Models (4/5 pages)

### Section 3.1 - Metropolis
* History of the method
* How it works
* Limitations
#### 3.1.1 - Critical Slowing Down
  * Explain this limitation better

### Section 3.2 - Random Path Sampling
* How the method works
* Its limitations

### Section 3.3 - Wang-Landau
* How the method works
* Its improvements relative to the Metropolis 
* Its problems 
* How some people tried solving them

## Plots for Chapter 3
* Schemes for the RPS mehtod (like in the JA ppt)
* Incomplete JDOS from the RPS method to show its weakness
* Plotted JDOS for the Wang-Landau sampling
* Fig.4 from the 1/t paper (2007)
* Fig.2 from the convergence paper (2005)

-----------------------------------------

## Chapter 4 - Introduction to the Flat Scan Sampling (5/7 pages)
* Introduce the FSS from the problems of the WL

### Section 4.1 - Background
* How the method was created
* Its relation to RPS and WL
#### 4.1.1 - Algorithm
  * How it actually works
### Section 4.2 - Implementation for the Ising model
#### 4.2.1 - Single Core
  * Single Core C++
  * Efforts to make the simulation more accurate -> skip
#### 4.2.2 - Performance
  * Comparison of Wall times for different lattices and sizes (for all implementations, MATLAB inc.)
#### 4.2.3 - MPI
  * One manager and N-1 workers (probably not included)
  * One manager and N workers
  * Efforts to make the simulation more accurate -> WL_shuffle
#### 4.2.4 - Multicore Scaling
  * Talk about Amdhal's Law and how it fits to the results
  * Multi-Core scalling of the method

## Plots for Chapter 4
* Scheme of how it works
* Compare computed magnetisation with one from the WL and Metropolis for SC
* Results from skip and WL_shuffle
* Table from single core benchmarks
* Plot of the theoretical Amdhal's Law
* Plot of scaling
* Wall time plots with REP 

-----------------------------------------

## Chapter 5 - Validation and Convergence of the FSS (4/6 pages)
* Plots of statistical analysis and REP analysis

### Section 5.1 - Validation
* Compare thermodynamic results with Metropolis 
  * <|M|> (met) vs MmF (fss) for large L (met) and reasonalbe L (fss)
  * <|M|> for metropolis and fss same L
* Compare Tc values
* Compare computed JDOS with exact JDOS L4
* Compare DOS with Beale solution L16 (?) (just sum the magnetisations)

### Section 5.2 - Convergence
* Statistical analysis of REP L8
* Error for the L4 with the exact JDOS

## Plots for Chapter 5
* 2 magnetisation plots
* 1 plot with the TC from different lattices
* plto that compares the Beale exact DOS (?)
* plot that compares the estimated JDOS with the L4 exact solution
* 2 plots from the slice and checkerboard configs for the JDOS L8
* 1/2 plots of the deviation from the exact JDOS L4
* var vs 1/sqrt(REP) plot (linear relation)

------------------------------------------

## Chapter 6 - Wang-Landau Comparison (3/5 pages)
* Why compare with wang landau?
* Briefly explaing the strategy of both methods
* Advantages and disadvantages from both methods

### Section 6.1 - Statistical Analysis
* Error plots for the different configurations L8 or L16 with mc results
* Saturation of the error even for large values of flatness
* Should I reproduce the results from the 1/t method to compare? Or introduce the S parameter? 

### Section 6.2 - Estimating Tc for infinite lattice
* According to the 1/t paper the WL Tc(inf) is wrong
* A comparison of this should be good

## Plots for Chapter 6
* 2/3 plots to compare config mean and var 
* 1 plot L4 exact JDOS comparison
* 2 Tc extrapolation plot with the fss and wl results ss, sc
* wall time plots (WL is random)

------------------------------------------

## Chapter 7 - Finite Size Scalling (3/- pages)


------------------------------------------

## Chapter 8 - Conclusion and Future Work (1 page)
------------------------------------------

And if there is time left talk about the SpinS Ising Model

