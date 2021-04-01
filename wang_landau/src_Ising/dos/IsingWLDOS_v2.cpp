/*
    Wang Landau Sampling for the 2D Ising Model
    João Inácio, Sep. 1, 2020
*/
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "VectorFunctions.h"
#include "IsingFunctions.h"

void normalizeDOS(double*, int);

const int J = 1;                // Interaction strength
const int L = 2;                // Lattice size
const int N_SPINS = L * L;      // Number of particles

int main() {
    // Start measuring time
    clock_t start = clock();

    double* aux = createVector(0, L-1, 1);
    double* minus = shift(aux, L, -1);
    double* plus = shift(aux, L, 1);

    double f = exp(1);          // Modification factor
    double p = 0.90;            // Flatness -> min(Hist) > avg(Hist)*p

    // All of the possible configurations
    int NE = N_SPINS - 1;
    int energies[NE];
    generateEnergies(energies, N_SPINS);
    
    // Create spins
    int spins[L * L];
    randSpins(spins, L);
    // printSpins(spins, L);
    int E = CEnergy(spins, plus, minus, J, L, N_SPINS); // Energy of the first configuration
    auto itr = std::find(energies, energies + NE, E);
    int idxE = std::distance(energies, itr);

    // DOS array and Histogram
    double lngE[NE];
    setZero(lngE, NE);
    double hist[NE];
    setZero(hist, NE);
    
    // MCSweeps counter
    int MCSweeps = 0;

    // Random seed
    // srand((int)time(0));

    while(f > (1 + pow(10, -8))) {
        // Each loop is one MCSweep
        for (int ni = 0; ni < N_SPINS; ni++) {
            int i = (rand() % L); // Generates a random int (1 - L)
            int j = (rand() % L);

            // dE = 2 * J * S * sum(neighbours) and dE = ENew - E, so
            int S = spins[i * L + j];
            int sumN = sumNeigh(i, j, spins, plus, minus, L, N_SPINS);
            int ENew = E + 2 * J * S * sumN;

            auto itr = std::find(energies, energies + NE, ENew);
            int idxENew = std::distance(energies, itr);

            // Flip the spin
            double ratio = exp(lngE[idxE] - lngE[idxENew]);
            double P = std::min(ratio, 1.0);

            if (P > ((double) rand() / (RAND_MAX))) {
                spins[i * L + j] = - S;
                E = ENew;
                idxE = idxENew;
            }
            // printSpins(spins);
            // Update the histogram and g(E)
            hist[idxE] += 1;
            lngE[idxE] += log(f);
            // print1D(lngE, NE);
        }

        MCSweeps += 1;

        // Check at each 10000 sweeps
        if (MCSweeps % 10000 == 0) {
            double avgH = average(hist, NE);
            int minH = *std::min_element(hist, hist + NE);

            for (int i = 0; i < NE; i++) 
                std::cout << energies[i] << " : " << hist[i] << std::endl;

            if (minH > avgH * p) {
                printf("%d: the histogram is flat. Min: %d Avg: %.1f f: %.8f\n", MCSweeps, minH, avgH, f);
                setZero(hist, NE);
                f = sqrt(f);
                MCSweeps = 0;
            }
        }
    }

    // Normalize the DOS
    normalizeDOS(lngE, NE);

    // Get the actual DOS
    double g[NE];

    for (int i = 0; i < NE; i++) {
        g[i] = exp(lngE[i]);
    }

    print1D(g, NE);

    // Writing the DOS in a file
    std::string fileName = "DOS_" + std::to_string(L) + "L_WL_CPP.txt";

    std::ofstream file;
    file.open(fileName);
    
    for (int i = 0; i < NE; i++) {
        file << lngE[i] << " ";
    }

    file.close();

    // Erase arrays from the heap
    delete[] aux, minus, plus;

    // Stop measuring time and calculate the elapsed time
    clock_t end = clock();
    double elapsed = double(end - start)/CLOCKS_PER_SEC;
    printf("Time measured: %.7f seconds.\n", elapsed);
    return 0;
}

void normalizeDOS(double* lngE, int NE) {
    double lngE0 = lngE[0];
    for (int i = 0; i < NE; i++)  
        lngE[i] = lngE[i] - lngE0 + log(2);
}
