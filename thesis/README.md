# Report Outline (30 pages)
---

# NEW OUTLINE

---

## Chapter 1 - Ferromagnetism and the Ising Model

### 1.1 - Ferromagnetism
* What is ferromagnetism
* What materials display it
* What does it have to do with computational physics
* Good models that model it

### 1.2 - Ising Model
* What is the Ising model
* History
* Phase transitions
#### 1.2.1 - Joint Density of States
* What is the DOS and JDOS
* DOS vs JDOS
#### 1.2.2 - Thermodynamics
* Thermodynamic relations that can be computed from the DOS and JDOS
#### 1.2.3 - Relevance
* Where can we find the model

## Plots for Chapter 1
* :heavy_check_mark: Ferromagnet/Paramagnet scheme 
* :heavy_check_mark: Ising phase transition plot/scheme
* :heavy_check_mark: JDOS table for a L2_SS system
* :heavy_check_mark: Plotted JDOS_exact for L4_SS
---

## Chapter 2 - Monte Carlo Methods Applied to the Ising Model
* Why solve this models by computer simulations
* Why are stochastic methods good to simulate these models

### 2.1 - Metropolis
* History and its legacy
* How it works
* Limitations
#### 2.1.1 - Critical Slowing Down
 * Explain this limitation a little better

### 2.2 - Wang-Landau
* Efforts in the 90s to find better methods than Metropolis
#### 2.2.1 - Algorithm
* How it works
#### 2.2.2 - Variations
* 1/t
* S parameter
#### 2.2.3 - Success and Limitations
* How it solved some problems and brought new ones

* Applications

## Plots for Chapter 2
* Something for the Metropolis
* DOS and JDOS computed by the WL samping
* Error saturation with standart WL as a function of f_final and flatness
* var vs something to show that is not linear
* Fig.4 from the 1/t paper (2007)
* Fig.2 from the convergence paper (2005)
---

## Chapter 3 - Flat Scan Sampling
* Intoduce FSS from the problems from the WL method

### 3.1 - Background 
* Explain RPS method 
* Relation to both the RPS and WL

### 3.2 - Algorithm
* How it works

### 3.3 - Implementations
#### 3.3.1 - Single Core
 * Single core C++
 * Explain the thought process
 * Efforts to make the simulation better -> skip
#### 3.3.2 - MPI
 * N worker implementation
 * Explain implementation
 * Efforts to make the simulation better -> WL shuffle

### 3.4 - Validation and Convergence
* Validation for L4_SS
* rep analysis for convergence
* Show that the method converges to the exact solution with and error proportional to 1/sqrt(rep)
* Show the same for L8 and L16 SS systems

### 3.5 - Performance
* Single core performance -> linear with REP
* Multi core performance 
#### 3.5.1 - Amdhal's Law and Parallel Scaling
 * Amdhal's law 
 * Fit to results
 * Estimation of speedup for infinite cores

### 3.6 - Comparison with Wang-Landau Sampling
* Why compare with wang landau?
* Briefly explaing the strategy of both methods
* Advantages and disadvantages from both methods
#### 3.6.1 - Statistical Analysis
 * Error plots for the different configurations L8 or L16 with mc results
 * Saturation of the error even for large values of flatness
 * Reproduce the results for 1/t and S parameter variations to compare

## Plots for Chapter 3
* Scheme for the RPS
* JDOS to show that the RPS sometimes doesn't find all of the pairs (E, M)

### Implementations
* Plot to show the effects of skip
* Scheme for the MPI implmentation (?)
* Plot to show the effects of WL shuffle

### Validation and Convergence
* :heavy_check_mark: rep analysis plots to show that the error becomes 0 and converges to the right solution when we increase rep
* :heavy_check_mark: var vs 1/sqrt(rep) to show that the method converges to the exact solution while being inear with 1/sqrt(rep)
* Show the same for L8 and L16 systems
* :heavy_check_mark: Plot that show the exact and computed DOS for q=9 (L4) and the error in each point
* :heavy_check_mark: Show the histogram fit for the mean error in the L4 system
* :heavy_check_mark: Show the histogram fit for the 3 configs in the L4 system

### Performance
* :heavy_check_mark: Plot of the wall time vs rep
* :heavy_check_mark: Same plot for the q time per E
* Plot of the wall time vs rep for MPI
* :heavy_check_mark: Theoretical behaviour for the Amdhal's law
* :heavy_check_mark: Plot for the 1/n to extrapolate the S(inf)
* :heavy_check_mark: Plot for Amdhal fit and S(inf) with the data points

### Comparison with Wang-Landau Sampling
* Error comparison with WL as a function o f_final and flatness
* Compare error with the other implementations
* Wall time comparison
* Try to reproduce the plots from the S and 1/t papers

---

## Chapter 4 - Thermodynamics and Finite Size Scaling



## Plots for Chapter 4
---

## Chapter 5 - Conclusion and Future Work
* Give a little peak at spinS computations

## Plots for Chapter 5
* Maybe a plot of something from spinS Ising
---

# OLD OUTLINE

---

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

---

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

---

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

---

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

---

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

---

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

---

## Chapter 7 - Finite Size Scalling (3/- pages)


---

## Chapter 8 - Conclusion and Future Work (1 page)
---

And if there is time left talk about the SpinS Ising Model

