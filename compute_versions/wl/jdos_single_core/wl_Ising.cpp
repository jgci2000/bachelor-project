//
// Wang Landau sampling for the Ising 1/2 Model 
// João Inácio, Mar. 30th, 2021
//
// This version is single core and makes use the Ising class
//

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdint.h>
#include <chrono>
#include <cmath>
#include <array>
#include <vector>
#include <map>

#include "Ising.h"
#include "WL_Functions.h"
#include "splimix64.h"
#include "xoshiro256++.h"

using std::cout;
using std::endl;
using std::vector;
using std::array;
using std::string;

#define ll          long long
#define ld          long double


// Seed for the RNG. If SEED equal 0, use random seed.
#define SEED        0

// Size of the Ising Lattice
#define L_LATTICE   8
// LATTICE_NUM -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex 
#define LATTICE_NUM 1

// Output location
#define SAVE_DIR    "./data/"

int main(int argc, char **argv)
{
    // Set the seed for xoshiro256++
    uint64_t seed = SEED;
    if (seed == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }

    splitmix64_seed(seed);
    for (int i = 0; i < 4; i++)
        s[i] = splitmix64();
    
    // Initialize Ising and set parameters for FSS computations
    Ising ising(L_LATTICE, LATTICE_NUM);

    int q_max = (ising.NM + 1) / 2 - 2;

    if (ising.NM % 2 == 0)
        q_max = ising.NM / 2 - 3;

    int run = atoi(argv[1]);
    double f = exp(1);
    double f_final = 1 + pow(10, - atoi(argv[2]));
    double flatness = 0.90;

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_" + std::to_string(ising.NN) + "NN_L" + std::to_string(ising.L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + std::to_string(ising.N_atm) + "d2.txt";
    string save_file = std::to_string(run) + "_JDOS_WL_Ising_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_L" + std::to_string(ising.L) + "_f" + std::to_string((int) - log10(f_final - 1)) + "_flatness" + std::to_string((int) (flatness * 100));

    ising.read_NN_talbe(NN_table_file_name);
    ising.read_norm_factor(norm_factor_file_name);

    ld *ln_JDOS = new ld[ising.NE * ising.NM];
    ld *JDOS = new ld[ising.NE * ising.NM];
    ll *hist = new ll[ising.NE * ising.NM];

    for (int i = 0; i < ising.NE * ising.NM; i++)
    {
        JDOS[i] = 0;
        ln_JDOS[i] = 0;
        hist[i] = 0;
    }

    ll mc_sweep = 0;

    // Random spins config
    for (int i = 0; i < ising.N_atm; i++) 
    {
        if ((rand_xoshiro256pp() % 2) + 1 == 1) 
            ising.spins_vector[i] = + 1;
        else ising.spins_vector[i] = - 1;
    }

    int E_config, M_config = 0;
    for (int i = 0; i < ising.N_atm; i++)
    {
        for (int a = 0; a < ising.NN; a++)
            E_config += - ising.spins_vector[i] * ising.spins_vector[ising.NN_table[i * ising.NN + a]];

        M_config += ising.spins_vector[i];
    }
    E_config = E_config / 2;

    ising.set_E_config(E_config);
    ising.set_M_config(M_config);

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "L: " + std::to_string(ising.L) + " | f_final: 1+1E" + std::to_string((int) log10(f_final - 1)) + " | flatness: " + std::to_string((int) (flatness * 100)) + " | dim: " + std::to_string(ising.dim) + "D | lattie: " + ising.lattice;
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // Wang Landau Sampling
    std::chrono::_V2::steady_clock::time_point loop_start;

    // Main loop
    while(f > f_final)
    {
        if (mc_sweep == 0)
            loop_start = std::chrono::steady_clock::now();
        
        for (int idx = 0; idx < ising.N_atm; idx++) 
        {
            int flip_idx = (rand_xoshiro256pp() % ising.N_atm);

            int delta_E = 0;
            for (int a = 0; a < ising.NN; a++)
                delta_E += - ising.spins_vector[flip_idx] * ising.spins_vector[ising.NN_table[flip_idx * ising.NN + a]];
            
            int new_E_config = ising.E_config - 2 * delta_E;
            int new_M_config  = ising.M_config - 2 * ising.spins_vector[flip_idx];
            int new_idx_E_config = ising.energies[new_E_config];
            int new_idx_M_config = ising.magnetizations[new_M_config];
            
            ld ratio = exp(ln_JDOS[ising.idx_E_config * ising.NM + ising.idx_M_config] - ln_JDOS[new_idx_E_config * ising.NM + new_idx_M_config]);

            if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio)
            {
                ising.spins_vector[flip_idx] = - ising.spins_vector[flip_idx];
                ising.set_E_config(new_E_config);
                ising.set_M_config(new_M_config);
            }

            hist[ising.idx_E_config * ising.NM + ising.idx_M_config]++;
            ln_JDOS[ising.idx_E_config * ising.NM + ising.idx_M_config] += log(f);
        }

        mc_sweep++;

        if (mc_sweep % 10000 == 0) 
        {
            long double avg_h = average_hist(hist, ising.NE * ising.NM);
            int min_h = min_hist(hist, ising.NE * ising.NM);
            
            if (min_h >= avg_h * flatness)
            {
                auto loop_end = std::chrono::steady_clock::now();
                double loop_dur = (double) (std::chrono::duration_cast<std::chrono::microseconds> (loop_end - loop_start).count()) * pow(10, -6);

                now = time(0);
                t = ctime(&now); t.pop_back();
                
                string console_output = t + " | f: 1+1E" + std::to_string(log10(f - 1)) + "/1+1E" + std::to_string((int) log10(f_final - 1)) + " | sweeps: " + std::to_string(mc_sweep) + " | flat time: " + std::to_string(loop_dur) + "s";
                string data_line = "1+1E" + std::to_string(log10(f - 1)) + " 1+1E" + std::to_string((int) log10(f_final - 1)) + " " + std::to_string(loop_dur) + " " + std::to_string(mc_sweep) + " " + std::to_string(min_h)
                + " " + std::to_string(avg_h);

                console_log.push_back(console_output);
                data.push_back(data_line);

                cout << console_output << endl;

                f = sqrt(f);
                mc_sweep = 0;

                for (int i = 0; i < ising.NE * ising.NM; i++)
                    hist[i] = 0;
            }
        }
    }

    // Normalize JDOS
    for (int i = 0; i < ising.NE; i++)
    {
        for (int j = 0; j < ising.NM; j++)
        {
            if (ln_JDOS[i * ising.NM + j] > 0) 
                JDOS[i * ising.NM + j] = exp(ln_JDOS[i * ising.NM + j] - ln_JDOS[0] + log(2)) / 2; 
        }
    }

    // Stop mesuring time
    auto method_end = std::chrono::steady_clock::now();
    double runtime = (double) (std::chrono::duration_cast<std::chrono::microseconds> (method_end - method_start).count()) * pow(10, -6);
    now = time(0);
    t = ctime(&now); t.pop_back();

    cout << endl;
    cout << "Runtime: " << std::setw(8) << runtime << " seconds." << endl;
    cout << "Simulation ended at: " << t << endl;

    // Write JDOS to file
    std::ofstream file1((string) SAVE_DIR + save_file + ".txt");
    for (int i = 0; i < ising.NE; i++) 
    {
        for (int j = 0; j < ising.NM; j++) 
            file1 << JDOS[i * ising.NM + j] << " ";
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2((string) SAVE_DIR + save_file + "_data.txt");
    file2 << "f f_max loop_dur mc_sweeps min_h avg_h \n"; 
    for (int i = 0; i < data.size(); i++)
        file2 << data[i] << "\n";
    file2 << runtime << "\n";
    file2.close();

    std::ofstream file3((string) SAVE_DIR + save_file + "_console_logs.txt");
    for (int i = 0; i < console_log.size(); i++)
        file3 << console_log.at(i) << "\n";
    file3.close();

    delete[] JDOS, ln_JDOS, hist;

    return 0;
}
