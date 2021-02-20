/*
    New Random Path Sampling for the 2D Ising Model with Parallelization
    João Inácio, Nov. 2, 2020
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdint.h>
#include <chrono>
#include <omp.h>


#include "./.includes/VectorFunctions.h"
#include "./.includes/IsingFunctions.h"


using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
//
//  New Random Path Sampling for the 2D Ising Model using openMP - v4
//  João Inácio, Dec. 18, 2020
//
//


#define WRITE_FILE                  1           // Enables writing the JDOS to a text file
#define SEED                        0
#define NUMBER_THREADS_REQUESTED    1           // Select the number of theards to request. The actual number of threads use might be lower.
#define MAX_THREADS                 omp_get_max_threads()
#define lli                         long long int


const int J = 1;                    // Interaction strength
const int L = 4;                    // Lattice size
const int N_SPINS = L * L;          // Number of particles
const int NN = 4;                   // Number of neighbours

const int q_max = N_SPINS / 2 + 1;  // Index of magnetization and the number os spins down
const int skip = N_SPINS;           // Scan a configuration each skip times to get more random results
const lli REP = (int) pow(10, 3);   // Number of sampled configurations in one energy

const int max_E = (1.0 / 2.0) * NN * N_SPINS;
const int max_M = N_SPINS;

const int NE = 1 + (max_E / 2);     // Number of allowed energies
const int NM = N_SPINS + 1;         // Number of allowed magnetizations


// Some functions
struct xorshift64s_state {
    uint64_t a;
};

uint64_t xorshift64s(struct xorshift64s_state *state)
{
    uint64_t x = state->a;	
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    state->a = x;
    return x * UINT64_C(0x2545F4914F6CDD1D);
}


int main(int argc, char **argv) 
{
    // Set a seed for the random number generator
    uint64_t seed;
    if (SEED == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }
    else
        seed = SEED;
    
    xorshift64s_state state = {.a = seed};
    
    // Auxiliary vectors and lattice
    array<int, N_SPINS> lattice;
    create_vector(lattice, 0, N_SPINS - 1, 1);
    
    array<int, N_SPINS> nnxpos = circshift(lattice, 0, - 1);
    array<int, N_SPINS> nnxneg = circshift(lattice, 0, 1);
    array<int, N_SPINS> nnypos = circshift(lattice, 1, 0);
    array<int, N_SPINS> nnyneg = circshift(lattice, - 1, 0);

    array<int, NE> energies;
    array<int, NM> magnetizations;
    create_vector(energies, - max_E, max_E, 4);
    create_vector(magnetizations, - max_M, max_M, 2);

    double *JDOS = new double[NE * NM];
    double *JDOS_frac = new double[NE * NM];

    for (int i = 0; i < NE * NM; i++)
    {
        JDOS[i] = 0;
        JDOS_frac[i] = 0;
    }

    JDOS[0] = 1;
    JDOS[NM - 1] = 1;
    JDOS_frac[0] = 1;
    JDOS_frac[NM - 1] = 1;

    lli *neo_previous = new lli[NE * NE];
    
    // Open the normalization factor
    std::ifstream norm_factor_file("./.norm_factors/norm_factor_Ising_Natm_" + std::to_string(N_SPINS) + ".txt");
    string line;
    long double *norm_factor = new long double[NM];

    if (norm_factor_file.is_open()) 
    {
        for (int i = 0; std::getline(norm_factor_file, line); i++)
        {
            std::istringstream iss(line);
            iss >> norm_factor[i];
        }
        norm_factor_file.close();
    }
    else 
    {
        cout << "Unable to open file. Invalid lattice size. Or the file isn't on the correct directory. Please read README.txt" << endl;
        return 0;
    }

    // Defining the threads number
    if (NUMBER_THREADS_REQUESTED >= MAX_THREADS) omp_set_num_threads(MAX_THREADS);
    else omp_set_num_threads(NUMBER_THREADS_REQUESTED);

    int number_working_threads = omp_get_max_threads();
    int REP_thread = REP / number_working_threads;
    
    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();
    
    string console_output = "L: " + std::to_string(L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | threads: " + std::to_string(number_working_threads) + " | REP per thread " + std::to_string(REP_thread);
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // New RPS algorithm
    for (int q = 0; q < q_max - 1; q++) 
    {
        auto q_start = std::chrono::steady_clock::now();
        
        // Global variables
        for (int i = 0; i < NE * NE; i++)
            neo_previous[i] = 0;

        lli accept_counter = 1;
        lli reject_counter = 0;
        lli k = 0;
        lli k_saved = 0;
        lli sample_ratio = 0;

        // Parallelization        
        #pragma omp parallel //shared(neo_previous, accept_counter, reject_counter, k, k_saved, sample_ratio)
        {
            int ID = omp_get_thread_num();

            array<double, NE> WL_log_DOS;
            for (int i = 0; i < NM; i++)
                WL_log_DOS[i] = log(JDOS[i * NM + q]);

            array<int, NE> hist_WL_thread; hist_WL_thread.fill(0);
            array<int, NE> hist_E_selected_thread; hist_E_selected_thread.fill(0);

            array<int, N_SPINS> spins_WL; spins_WL.fill(1);
            int E_WL_old = - max_E;

            vector<int> pos;
            vector<int> neg;
            for (int i = 0; i < lattice.size(); i++)
                pos.push_back(lattice[i]);

            if (q >= 1)
            {
                for (int idx = 0; idx < q; idx++)
                {
                    int idx_pos = xorshift64s(&state) % pos.size();
                    int flipped_pos = pos.at(idx_pos);

                    pos.erase(pos.begin() + idx_pos);
                    neg.push_back(flipped_pos);

                    spins_WL[flipped_pos] = - 1;

                    int delta_E = - J * spins_WL[flipped_pos] * 
                    (spins_WL[nnxpos[flipped_pos]] + 
                    spins_WL[nnxneg[flipped_pos]] + 
                    spins_WL[nnypos[flipped_pos]] + 
                    spins_WL[nnyneg[flipped_pos]]);
                    E_WL_old += 2 * delta_E;
                }
            }

            int idx_E_WL_old = binary_search(energies, E_WL_old);
            int E_old = E_WL_old;
            int idx_E_old = idx_E_WL_old;

            for (int i = 0; i < pos.size(); i++)
            {
                int flipped_pos_scan = pos.at(i);
                spins_WL[flipped_pos_scan] = - 1;

                int delta_E = - J * spins_WL[flipped_pos_scan] * 
                    (spins_WL[nnxpos[flipped_pos_scan]] + 
                    spins_WL[nnxneg[flipped_pos_scan]] + 
                    spins_WL[nnypos[flipped_pos_scan]] + 
                    spins_WL[nnyneg[flipped_pos_scan]]);
                int E_new = E_old + 2 * delta_E;

                int idx_E_new = binary_search(energies, E_new);

                spins_WL[flipped_pos_scan] = 1;

                #pragma omp atomic update
                    neo_previous[idx_E_old * NE + idx_E_new]++;
            }

            hist_WL_thread[idx_E_WL_old]++;
            hist_E_selected_thread[idx_E_WL_old]++;

            lli k_saved_thread = 1;
            lli k_thread = 2;
            int c = 0;
            while (min_hist(hist_E_selected_thread) < REP_thread)
            {
                int idx_pos = xorshift64s(&state) % pos.size();
                int flipped_pos = pos.at(idx_pos);
                pos.erase(pos.begin() + idx_pos);
                neg.push_back(flipped_pos);

                spins_WL[flipped_pos] = - 1;

                int delta_E = - J * spins_WL[flipped_pos] * 
                    (spins_WL[nnxpos[flipped_pos]] + 
                    spins_WL[nnxneg[flipped_pos]] + 
                    spins_WL[nnypos[flipped_pos]] + 
                    spins_WL[nnyneg[flipped_pos]]);
                int E_WL_new = E_WL_old + 2 * delta_E;

                int idx_neg = xorshift64s(&state) % neg.size();
                int flipped_neg = neg.at(idx_neg);
                neg.erase(neg.begin() + idx_neg);
                pos.push_back(flipped_neg);

                spins_WL[flipped_neg] = 1;

                delta_E = - J * spins_WL[flipped_neg] * 
                    (spins_WL[nnxpos[flipped_neg]] + 
                    spins_WL[nnxneg[flipped_neg]] + 
                    spins_WL[nnypos[flipped_neg]] + 
                    spins_WL[nnyneg[flipped_neg]]);
                E_WL_new += 2 * delta_E;
                int idx_E_WL_new = binary_search(energies, E_WL_new);

                double ratio = exp(WL_log_DOS[idx_E_WL_old] - WL_log_DOS[idx_E_WL_new]);

                if (ratio >= 1 || (ratio * 10000) > ((double) (xorshift64s(&state) % 10000)))
                {
                    E_WL_old = E_WL_new;
                    idx_E_WL_old = idx_E_WL_new;
                    hist_WL_thread[idx_E_WL_old]++;

                    #pragma omp atomic update
                        accept_counter++;

                    if (k_thread >= k_saved_thread + skip && hist_E_selected_thread[idx_E_WL_old] <= REP_thread) 
                    {
                        hist_E_selected_thread[idx_E_WL_old]++;
                        k_saved_thread = k_thread;

                        int E_old = E_WL_old;
                        int idx_E_old = idx_E_WL_old;

                        for (int i = 0; i < pos.size(); i++)
                        {
                            int flipped_pos_scan = pos.at(i);
                            spins_WL[flipped_pos_scan] = - 1;

                            int delta_E = - J * spins_WL[flipped_pos_scan] * 
                                (spins_WL[nnxpos[flipped_pos_scan]] + 
                                spins_WL[nnxneg[flipped_pos_scan]] + 
                                spins_WL[nnypos[flipped_pos_scan]] + 
                                spins_WL[nnyneg[flipped_pos_scan]]);
                            int E_new = E_old + 2 * delta_E;

                            int idx_E_new = binary_search(energies, E_new);

                            spins_WL[flipped_pos_scan] = 1;

                            #pragma omp atomic update
                                neo_previous[idx_E_old * NE + idx_E_new]++;
                        }
                    }
                }
                else
                {   
                    spins_WL[flipped_neg] = - 1;
                    spins_WL[flipped_pos] = 1;
                    hist_WL_thread[idx_E_WL_old]++;

                    #pragma omp atomic update
                        reject_counter++;

                    if (k_thread >= k_saved_thread + skip && hist_E_selected_thread[idx_E_WL_old] <= REP_thread) 
                    {
                        hist_E_selected_thread[idx_E_WL_old]++;
                        k_saved_thread = k_thread;

                        int E_old = E_WL_old;
                        int idx_E_old = idx_E_WL_old;

                        for (int i = 0; i < pos.size(); i++)
                        {
                            int flipped_pos_scan = pos.at(i);
                            spins_WL[flipped_pos_scan] = - 1;

                            int delta_E = - J * spins_WL[flipped_pos_scan] * 
                                (spins_WL[nnxpos[flipped_pos_scan]] + 
                                spins_WL[nnxneg[flipped_pos_scan]] + 
                                spins_WL[nnypos[flipped_pos_scan]] + 
                                spins_WL[nnyneg[flipped_pos_scan]]);
                            int E_new = E_old + 2 * delta_E;

                            int idx_E_new = binary_search(energies, E_new);

                            spins_WL[flipped_pos_scan] = 1;

                            #pragma omp atomic update
                                neo_previous[idx_E_old * NE + idx_E_new]++;
                        }
                    }
                }
                k_thread++;

                // for (int i = 0; i < NE; i++)
                //     cout << hist_E_selected_thread[i] << " ";
                // cout << endl;

                // for (int i = 0; i < NE; i++)
                // {
                //     for (int j = 0; j < NM; j++)
                //         cout << neo_previous[i * NM + j] << " ";
                //     cout << endl;
                // }
                if (q == 2 && (c == 10000000 || c == 0))
                {
                    for (int i = 0; i < pos.size(); i++)
                        cout << pos.at(i) << " ";
                    cout << endl;

                    for (int i = 0; i < neg.size(); i++)
                        cout << neg.at(i) << " ";
                    cout << endl;

                    for (int i = 0; i < N_SPINS; i++)
                        cout << spins_WL[i] << " ";
                    cout << endl;
                }
                c++;
            }

            int hits = 0;
            for (int i = 0; i < NE; i++)
                if (JDOS[i * NM + q] > 0)
                    hits++;
            
            lli sample_ratio_thread = REP_thread * hits * skip / k_thread;
               
            #pragma omp atomic update
                k += k_thread;
            #pragma omp atomic update
                k_saved += k_saved_thread;
            #pragma omp atomic update
                sample_ratio += sample_ratio_thread;

            for (int i = 0; i < NE; i++)
                cout << hist_E_selected_thread[i] << " ";
            cout << endl;
        }
        
        int hits = 0;
        for (int i = 0; i < NE; i++)
            if (JDOS_frac[i * NM + q] > 0)
                hits++;

        for (int i = 0; i < NE * NE; i++)
        {
            if (neo_previous[i] != 0)
            {
                int idx_old = i / NE;
                int idx_new = i % NE;

                double sum_neo_previous_old = 0;
                for (int j = 0; j < NE; j++)
                    sum_neo_previous_old += neo_previous[idx_old * NE + j];

                JDOS_frac[idx_new * NM + (q + 1)] += JDOS_frac[idx_old * NM + q] * neo_previous[idx_old * NE + idx_new] / sum_neo_previous_old;
            }
        }

        for (int i = 0; i < NE; i++)
            JDOS[i * NM + (q + 1)] = JDOS_frac[i * NM + (q + 1)] * norm_factor[q + 1];
        
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < NM; j++)
                cout << neo_previous[i * NM + j] << " ";
            cout << endl;
        }
        
        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();
        
        string console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max - 2) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits) + "s";
        string data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        + " " + std::to_string(k) + " " + std::to_string(accept_counter) + " " + std::to_string(reject_counter) + " " + std::to_string(sample_ratio / number_working_threads);
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
        if (WRITE_FILE == 1) 
        {
            std::ofstream file1("./Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP.txt");
            for (int i = 0; i < NE; i++) 
            {
                for (int j = 0; j < NM; j++) 
                    file1 << JDOS[i * NM + j] << " ";
                file1 << "\n";
            }
            file1.close();

            std::ofstream file2("./Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP_data.txt");
            file2 << "q q_max q_time hits q_time/hits k accept_counter reject_counter sampled_ratio/numer_working_threads\n"; 
            for (int i = 0; i < data.size(); i++)
                file2 << data[i] << "\n";
            file2 << runtime << "\n";
            file2.close();

            std::ofstream file3("./Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP_console_logs.txt");
            for (int i = 0; i < console_log.size(); i++)
                file3 << console_log.at(i) << "\n";
            file3.close();
        }

        delete[] JDOS, JDOS_frac, norm_factor, neo_previous;

    return 0;
}


