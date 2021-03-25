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

#include "Ising.h"
#include "Fss_Functions.h"
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
#define L_LATTICE   4
// LATTICE_NUM -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex 
#define LATTICE_NUM 1

// Output location and name
#define SAVE_DIR    "./Data/"


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

    int skip = ising.N_atm;
    ll REP = pow(10, 4);

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_" + std::to_string(ising.NN) + "NN_L" + std::to_string(ising.L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + std::to_string(ising.N_atm) + "d2.txt";
    string save_file = "JDOS_FSS_Ising_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_L" + std::to_string(ising.L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip);

    ising.read_NN_talbe(NN_table_file_name);
    ising.read_norm_factor(norm_factor_file_name);

    ld *JDOS = new ld[ising.NE * ising.NM];
    ll *hist = new ll[ising.NE];
    ll *hist_E_selected = new ll[ising.NE];

    for (int i = 0; i < ising.NE * ising.NM; i++)
        JDOS[i] = 0;
    JDOS[0] = 1;   

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "L: " + std::to_string(ising.L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | dim: " + std::to_string(ising.dim) + "D | lattie: " + ising.lattice;
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // Flat Scan Sampling

    // Scan and Compute JDOS at q = 1    
    for (int flip_idx = 0; flip_idx < ising.N_atm; flip_idx++)
    {
        int E_tmp1 = 0;
        for (int a = 0; a < ising.NN; a++)
            E_tmp1 += - 1;
        
        int E_tmp2 = 0;
        for (int a = 0; a < ising.NN; a++)
            E_tmp2 += 1;
        
        int E_tmp3 = - ising.max_E - E_tmp1 + E_tmp2;
        int idx_E_tmp3 = binary_search(ising.energies, E_tmp3);

        JDOS[idx_E_tmp3 * ising.NM + 1] += JDOS[0];
    }

    int sum_JDOS = 0;
    for (int i = 0; i < ising.NE; i++)
        if (JDOS[i * ising.NM + 1] > 0)
            sum_JDOS += JDOS[i * ising.NM + 1];

    for (int i = 0; i < ising.NE; i++)
        JDOS[i * ising.NM + 1] = JDOS[i * ising.NM + 1] * ising.norm_factor[1] / sum_JDOS;

    // Main Loop    
    for (int q = 1; q <= q_max; q++)
    {
        auto q_start = std::chrono::steady_clock::now();

        for (int i = 0; i < ising.NE; i++)
        {
            hist[i] = 0; 
            hist_E_selected[i] = 0;
        }

        // Random config at q
        for (int i = 0; i < ising.N_atm; i++)
            ising.spins_vector[i] = 1;
        ising.set_E_config(- ising.max_E);

        for (int idx = 1; idx <= q; idx++)
        {
            
        }


        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        // string console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max - 2) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits) + "s | shuffle time: " + std::to_string(shuffle_time) + "s";
        // string data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        // + " " + std::to_string(k) + " " + std::to_string(accept_counter) + " " + std::to_string(reject_counter);
        
        // console_log.push_back(console_output);
        // data.push_back(data_line);

        // cout << console_output << endl;
    }
    








    // Stop mesuring time
    auto method_end = std::chrono::steady_clock::now();
    double runtime = (double) (std::chrono::duration_cast<std::chrono::microseconds> (method_end - method_start).count()) * pow(10, -6);
    now = time(0);
    t = ctime(&now); t.pop_back();

    cout << endl;
    cout << "Runtime: " << std::setw(8) << runtime << " seconds." << endl;
    cout << "Simulation ended at: " << t << endl;
    
    delete[] JDOS;

    return 0;
}
