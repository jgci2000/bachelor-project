//
// MPI implementation of the Flat Scan Sampling for the Ising 1/2 Model
// João Inácio, Mar. 30th, 2021
//
// This version is parallelized and makes use the Ising class
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
#include <mpi.h>

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

// Output location
#define SAVE_DIR    "./data/"


int main(int argc, char **argv)
{
    // Root, rank and size
    int root = 0;
    int rank;
    int size;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set the seed for xoshiro256++
    uint64_t seed = SEED;
    if (seed == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }

    splitmix64_seed(seed * rank);
    for (int i = 0; i < 4; i++)
        s[i] = splitmix64();

    // Initialize Ising and set parameters for FSS computations
    Ising ising(L_LATTICE, LATTICE_NUM);

    int q_max = (ising.NM + 1) / 2 - 2;

    if (ising.NM % 2 == 0)
        q_max = ising.NM / 2 - 3;

    int skip = ising.N_atm;
    ll REP = pow(10, 3);
    ll REP_worker = REP / (size - 1);

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_" + std::to_string(ising.NN) + "NN_L" + std::to_string(ising.L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + std::to_string(ising.N_atm) + "d2.txt";
    string save_file = "JDOS_FSS_Ising_" + std::to_string(ising.dim) + "D_" + ising.lattice + "_L" + std::to_string(ising.L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip);

    ising.read_NN_talbe(NN_table_file_name);

    ld *JDOS_worker;
    ld *JDOS_root;
    ll *hist;
    ll *hist_E_selected;

    int hits_root;
    int hits_worker = 0;

    int q;

    // Flat Scan Sampling
    if (rank == root)
    {
        // Root
        ising.read_norm_factor(norm_factor_file_name);

        JDOS_worker = new ld[ising.NE * ising.NM];
        JDOS_root = new ld[ising.NE * ising.NM];

        for (int i = 0; i < ising.NE * ising.NM; i++)
        {
            JDOS_root[i] = 0;
            JDOS_worker[i] = 0;
        }
        JDOS_root[0] = 1;

        // Start measuring time
        vector<string> console_log;
        vector<string> data;
        time_t now;
        string t;

        now = time(0);
        t = ctime(&now); t.pop_back();

        string console_output = "L: " + std::to_string(ising.L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | dim: " + std::to_string(ising.dim) + "D | lattie: " + ising.lattice + " | walkers: " + std::to_string(size - 1) + " | REP/walker: " + std::to_string(REP_worker);
        console_log.push_back(console_output);

        cout << endl;
        cout << console_output << endl;
        cout << "Starting time: " << t << endl << endl;

        auto method_start = std::chrono::steady_clock::now();

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
            int idx_E_tmp3 = ising.energies[E_tmp3];

            JDOS_root[idx_E_tmp3 * ising.NM + 1] += JDOS_root[0];
        }

        int sum_JDOS = 0;
        for (int i = 0; i < ising.NE; i++)
            if (JDOS_root[i * ising.NM + 1] > 0)
                sum_JDOS += JDOS_root[i * ising.NM + 1];

        for (int i = 0; i < ising.NE; i++)
            JDOS_root[i * ising.NM + 1] = JDOS_root[i * ising.NM + 1] * ising.norm_factor[1] / sum_JDOS;

        console_output = t + " | q: " + std::to_string(0) + "/" + std::to_string(q_max);
        console_log.push_back(console_output);

        cout << console_output << endl;

        // Main loop
        for (q = 1; q <= q_max; q++)
        {
            auto q_start = std::chrono::steady_clock::now();

            MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(JDOS_root, ising.NE * ising.NM, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
            
            MPI_Reduce(JDOS_worker, JDOS_root, ising.NE * ising.NM, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            MPI_Reduce(&hits_worker, &hits_root, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

            for (int j = 0; j <= q; j++)
                for (int i = 0; i < ising.NE; i++)
                    JDOS_root[i * ising.NM + j] = JDOS_root[i * ising.NM + j] / (size - 1);
            
            hits_root = hits_root / (size - 1);

            ld sum_JDOS = 0;
            for (int i = 0; i < ising.NE; i++)
                if (JDOS_root[i * ising.NM + q + 1] > 0)
                    sum_JDOS += JDOS_root[i * ising.NM + q + 1];

            for (int i = 0; i < ising.NE; i++)
                JDOS_root[i * ising.NM + q + 1] = JDOS_root[i * ising.NM + q + 1] * ising.norm_factor[q + 1] / sum_JDOS;

            auto q_end = std::chrono::steady_clock::now();
            double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

            now = time(0);
            t = ctime(&now); t.pop_back();

            console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits_root) + " | q_time/E: " + std::to_string(q_time / hits_root) + "s";
            string data_line = std::to_string(q) + " " + std::to_string(q_max) + " " + std::to_string(q_time) + " " + std::to_string(hits_root) + " " + std::to_string(q_time / hits_root);

            console_log.push_back(console_output);
            data.push_back(data_line);

            cout << console_output << endl;
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
                file1 << JDOS_root[i * ising.NM + j] << " ";
            file1 << "\n";
        }
        file1.close();

        std::ofstream file2((string) SAVE_DIR + save_file + "_data.txt");
        file2 << "q q_max q_time hits q_time/hits \n";
        for (int i = 0; i < data.size(); i++)
            file2 << data[i] << "\n";
        file2 << runtime << "\n";
        file2.close();

        std::ofstream file3((string) SAVE_DIR + save_file + "_console_logs.txt");
        for (int i = 0; i < console_log.size(); i++)
            file3 << console_log.at(i) << "\n";
        file3.close();

        delete[] JDOS_root, JDOS_worker;
    }
    else
    {
        // Workers
        JDOS_worker = new ld[ising.NE * ising.NM];
        hist = new ll[ising.NE];
        hist_E_selected = new ll[ising.NE];
        
        do
        {
            MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(JDOS_worker, ising.NE * ising.NM, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);

            for (int i = 0; i < ising.NE; i++)
            {
                hist[i] = 0;
                hist_E_selected[i] = 0;
            }

            // Random config at q
            for (int i = 0; i < ising.N_atm; i++)
                ising.spins_vector[i] = 1;
            ising.set_E_config(- ising.max_E);

            array<vector<int>, 2> flip_list;
            for (int i = 0; i < ising.N_atm; i++)
                flip_list[0].push_back(i);

            int E_config = ising.E_config;
            for (int idx = 1; idx <= q; idx++)
            {
                int idx_tmp = rand_xoshiro256pp() % flip_list[0].size();
                int flipped_idx = flip_list[0].at(idx_tmp);
                ising.spins_vector[flipped_idx] = - 1;

                flip_list[1].push_back(flipped_idx);
                flip_list[0].erase(flip_list[0].begin() + idx_tmp);

                int delta_E = 0;
                for (int a = 0; a < ising.NN; a++)
                    delta_E += - ising.spins_vector[flipped_idx] * ising.spins_vector[ising.NN_table[flipped_idx * ising.NN + a]];

                E_config += 2 * delta_E;
            }
            ising.set_E_config(E_config);

            // Update Histograms
            hist[ising.idx_E_config]++;
            hist_E_selected[ising.idx_E_config]++;

            // Scan the first config
            for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
            {
                int delta_E = 0;
                for (int a = 0; a < ising.NN; a++)
                    delta_E += ising.spins_vector[flip_list[0].at(flip_idx)] * ising.spins_vector[ising.NN_table[flip_list[0].at(flip_idx) * ising.NN + a]];

                int E_tmp = ising.E_config + 2 * delta_E;
                int idx_E_tmp = ising.energies[E_tmp];

                JDOS_worker[idx_E_tmp * ising.NM + q + 1] += JDOS_worker[ising.idx_E_config * ising.NM + q] / REP_worker;
            }

            ll k = 1;
            int *new_spins_vector = new int[ising.N_atm];
            bool accepted = false;

            // Where the magic happens
            while (min_hist(hist_E_selected, ising.NE) < REP_worker)
            {
                // Get a new random condig at magnetization q
                if (!accepted)
                    for (int i = 0; i < ising.N_atm; i++)
                        new_spins_vector[i] =  ising.spins_vector[i];
                int new_E_config = 0;
                int new_idx_E_config = 0;

                // Flip a positive spin to a negative
                int idx_tmp1 = rand_xoshiro256pp() % flip_list[0].size();
                int flipped_idx1 = flip_list[0].at(idx_tmp1);
                new_spins_vector[flipped_idx1] = - 1;

                flip_list[1].push_back(flipped_idx1);
                flip_list[0].erase(flip_list[0].begin() + idx_tmp1);

                int delta_E = 0;
                for (int a = 0; a < ising.NN; a++)
                    delta_E += - new_spins_vector[flipped_idx1] * new_spins_vector[ising.NN_table[flipped_idx1 * ising.NN + a]];
                new_E_config = ising.E_config + 2 * delta_E;

                // Flip a negative spin to a positive
                int idx_tmp2 = rand_xoshiro256pp() % flip_list[1].size();
                int flipped_idx2 = flip_list[1].at(idx_tmp2);
                new_spins_vector[flipped_idx2] = 1;

                flip_list[0].push_back(flipped_idx2);
                flip_list[1].erase(flip_list[1].begin() + idx_tmp2);

                delta_E = 0;
                for (int a = 0; a < ising.NN; a++)
                    delta_E += - new_spins_vector[flipped_idx2] * new_spins_vector[ising.NN_table[flipped_idx2 * ising.NN + a]];
                new_E_config = new_E_config + 2 * delta_E;

                // Wang Landau criteria
                new_idx_E_config = ising.energies[new_E_config];
                ld ratio = JDOS_worker[ising.idx_E_config * ising.NM + q] / JDOS_worker[new_idx_E_config * ising.NM + q];

                if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio)
                {
                    for (int i = 0; i < ising.N_atm; i++)
                        ising.spins_vector[i] = new_spins_vector[i];

                    ising.set_E_config(new_E_config);
                    accepted = true;
                }
                else
                {
                    if (flipped_idx1 != flipped_idx2)
                    {
                        flip_list[0].pop_back();
                        flip_list[0].push_back(flipped_idx1);

                        flip_list[1].pop_back();
                        flip_list[1].push_back(flipped_idx2);
                    }
                    accepted = false;
                }

                hist[ising.idx_E_config]++;

                // Scan configuration
                if (hist_E_selected[ising.idx_E_config] < REP_worker && k % skip == 0)
                {
                    for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
                    {
                        int delta_E = 0;
                        for (int a = 0; a < ising.NN; a++)
                            delta_E += ising.spins_vector[flip_list[0].at(flip_idx)] * ising.spins_vector[ising.NN_table[flip_list[0].at(flip_idx) * ising.NN + a]];

                        int E_tmp = ising.E_config + 2 * delta_E;
                        int idx_E_tmp = ising.energies[E_tmp];

                        JDOS_worker[idx_E_tmp * ising.NM + q + 1] += JDOS_worker[ising.idx_E_config * ising.NM + q] / REP_worker;
                    }

                    hist_E_selected[ising.idx_E_config]++;
                }

                k++;
            }

            delete[] new_spins_vector;

            hits_worker = 0;
            for (int i = 0; i < ising.NE; i++)
                if (JDOS_worker[i * ising.NM + q] > 0)
                    hits_worker++;

            MPI_Reduce(JDOS_worker, JDOS_root, ising.NE * ising.NM, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            MPI_Reduce(&hits_worker, &hits_root, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
        } while (q != q_max);

        delete[] JDOS_worker, hist, hist_E_selected;
    }

    MPI_Finalize();
    return 0;
}

