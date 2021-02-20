/*
    Wang Landau Sampling for the 3D Ising Model JDOS
    João Inácio, Oct. 2, 2020
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <fstream>


#include "../../CPP Includes/VectorFunctions.h"
#include "../../CPP Includes/IsingFunctions.h"
// #include "../CPP Includes/IsingPrints.h"


#define WRITE_FILE 1            // Enables writing the JDOS to a text file

const int J = 1;                // Interaction strength
const int L = 8;                // Lattice size
const int N_SPINS = L * L * L;  // Number of particles
const int L2 = L * L;
const int NN = 6;               // Number of neighbours

const int MAX_E = (1.0 / 2.0) * NN * N_SPINS;
const int MAX_M = N_SPINS;

const int NE = 1 + (MAX_E / 2); // Number of allowed energies
const int NM = N_SPINS + 1;     // Number of allowed magnetizations

// (i, j, k) -> j + i * L + k * L ^ 2

int main() 
{
    // Free parameters
    double f = exp(1);        // Modification factor
    double f_final = 1 + pow(10, -8);
    double flatness = 0.90;   // Flatness

    // Some auxiliary vectors
    std::array<int, L> aux;
    create_vector(aux, 0, L - 1, 1);
    std::array<int, L> minus;
    shift(minus, aux, -1);
    std::array<int, L> plus;
    shift(plus, aux, 1);
    
    // All of the possible configurations
    std::array<int, NE> energies;
    std::array<int, NM> magnetizations;

    create_vector(energies, - MAX_E, MAX_E, 4);
    create_vector(magnetizations, - MAX_M, MAX_M, 2);
    
    // Fisrt configuration
    std::array<int, N_SPINS> spins;
    rand_spins(spins);

    int E = comp_energy3D(spins, plus, minus, J);
    int M = comp_magnetization3D(spins);

    int idx_E = binary_search(energies, E);
    int idx_M = binary_search(magnetizations, M);

    // DOS array and Histogram
    std::array<double, NE * NM> ln_JDOS;
    std::array<int, NE * NM> hist;

    ln_JDOS.fill(0);
    hist.fill(0);

    int mc_sweeps = 0;

    // Random seed
    srand((unsigned) time(NULL));

    // Start measuring time
    clock_t method_start = clock();
    clock_t loop_start;

    while (f > f_final) 
    {
        if (mc_sweeps == 0)
            loop_start = clock();
        
        for (int ni = 0; ni < N_SPINS; ni++) 
        {
            int i = (rand() % L); // Generates a random int (1 - L)
            int j = (rand() % L);
            int k = (rand() % L);

            int s = spins[j + i * L + k * L2];

            int sum_nei = spins[j + plus[i] * L + k * L2] + spins[j + minus[i] * L + k * L2] + spins[plus[j] + i * L + k * L2] +
                 + spins[minus[j] + i * L + k * L2] + spins[j + i * L + plus[k] * L2] + spins[j + i * L + minus[k] * L2];
            int new_E = E + 2 * J * s * sum_nei;
            int new_M  = M - 2 * s;

            int idx_new_E = binary_search(energies, new_E);
            int idx_new_M = binary_search(magnetizations, new_M);

            double ratio = exp(ln_JDOS[idx_E * NM + idx_M] - ln_JDOS[idx_new_E * NM + idx_new_M]);

            if ((ratio >= 1) || (ratio > ((double) rand() / (RAND_MAX)))) 
            {
                spins[j + i * L + k * L2] = - s;
                E = new_E;
                M = new_M;
                idx_E = idx_new_E;
                idx_M = idx_new_M;
            }

            hist[idx_E * NM + idx_M] += 1;
            ln_JDOS[idx_E * NM  + idx_M] += log(f);
        }

        mc_sweeps += 1;

        if (mc_sweeps % 10000 == 0) 
        {
            double avg_h = average_hist(hist);
            int min_h = min_hist(hist);

            int idx_min = std::distance(hist.begin(), std::find(hist.begin(), hist.end(), min_h));

            // for (int i = 0; i < NE / 2; i++) 
            // {
            //     for (int j = 0; j < NM / 2; j++) 
            //     {
            //         std::cout << hist[i * NM + j] << "  ";
            //     }
            //     std::cout << std::endl;
            // }

            // std::cout << min_h << "  " << hist[0] << "  " << avg_h << std::endl;
            // std::cout << NE << "  " << NM << std::endl;

            // return 0;

            if (min_h > avg_h * flatness) 
            {
                clock_t loop_end = clock();
                double loop_dur = (double) (loop_end - loop_start) / CLOCKS_PER_SEC;
                
                printf("%d: the histogram is flat. Time elapsed: %.3fs Min: %d Avg: %.1f f: %.8f\n", mc_sweeps, loop_dur, min_h, avg_h, f);

                f = sqrt(f);
                mc_sweeps = 0;
                hist.fill(0);
            }
        }
    }

    // Stop mesuring time
    clock_t method_end = clock();

    // Get the JDOS
    std::array<double, NE * NM> JDOS;
    JDOS.fill(0);
    comp_JDOS(ln_JDOS, JDOS);

    // Writing the DOS in a file
    if (WRITE_FILE == 1) 
    {
        std::ofstream file("JDOS_3D_" + std::to_string(L) + "L_WL_CPP.txt");

        for (int i = 0; i < NE; i++) 
        {
            for (int j = 0; j < NM; j++) 
                file << JDOS[i * NM + j] << " ";
            file << "\n";
        }
        
        file.close();
    }

    // Calculate the elapsed time
    double elapsed = (double) (method_end - method_start) / CLOCKS_PER_SEC;
    printf("Time measured: %.7f seconds.\n", elapsed);
    
    return 0;
}


