#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdint.h>
#include <chrono>
#include <array>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "fss_functions.h"

using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;

#define ll long long

//////////////////////////////////////////////////////////////////////////////////////
//
//  Flat Scan Sampling for the 2D Ising Model using MPI - v2
//  João Inácio, Jan. 19, 2021
//
//

/* Seed for the rng. If SEED == 0, then the seed is random. */
#define SEED 0

/* LATTICE -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex */
#define LATTICE 1
/* DIM -> 1 - 2D; 2 - 3D */
#define DIM 1
/* S-Spins particles */
#define S 1/2
/* Number of spin projections */
#define SZ 2 * (S + 1)

/* Lattice size */
#define L 4
/* Ineteraction strength */
#define J 1



/* Saving directory for the JDOS file */
#define SAVE_DIR(lattice, REP) "./Data/" + lattice + "/L" + std::to_string(L) + "/" + std::to_string((int) log10(REP)) + "/"

#if DIM == 1 && LATTICE == 1
    /* Number of particles */
#   define N_SPINS L * L
    /* Number os nearest neughbours */
#   define NN 4
#   define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#   define NEIGH_FILE "./neighbour_tables/neighbour_table_2D_SS_4NN_L" + std::to_string(L) + ".txt"
#   define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_2D_SS_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "SS";
#endif

#if DIM == 2 && LATTICE > 1 && LATTICE  < 7
#   if LATTICE == 2
        /* Number of particles */
#       define N_SPINS L * L * L
        /* Number os nearest neughbours */
#       define NN 6
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_SC_6NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_SC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "SC";
#   endif
#   if LATTICE == 3
        /* Number of particles */
#       define N_SPINS 2 * L * L * L
        /* Number os nearest neughbours */
#       define NN 8
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_BCC_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_BCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "BCC";
#   endif
#   if LATTICE == 4
        /* Number of particles */
#       define N_SPINS 4 * L * L * L
        /* Number os nearest neughbours */
#       define NN 12
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_FCC_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_FCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "FCC";
#   endif
#   if LATTICE == 5
        /* Number of particles */
#       define N_SPINS 2 * L * L * L
        /* Number os nearest neughbours */
#       define NN 12
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_HCP_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_HCP_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "HCP";
#   endif
#   if LATTICE == 6
        /* Number of particles */
#       define N_SPINS L * L * L
        /* Number os nearest neughbours */
#       define NN 8
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_Hex_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_Hex_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
const string lattice = "Hex";
#   endif

#endif

const int q_max = N_SPINS / 2 + 1;      // Max index of the magnetization computed (N_SPINS / 2 + 1 for half of the JDOS)
const int skip = N_SPINS;               // Scan a configuration each skip times to get more random results
// const ll REP = pow(10, 4);              // Number of sampled configurations in one energy

const int max_E = (1.0 / 2.0) * NN * N_SPINS;
const int max_M = N_SPINS;

const int NE = 1 + (max_E / 2);         // Number of allowed energies
const int NM = N_SPINS + 1;             // Number of allowed magnetizations

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
    int root = 0;

    array<int, NE> energies;
    array<int, NM> magnetizations;

    double *JDOS;
    double *JDOS_frac;
    int q;

    ll accept_counter = 0;
    ll reject_counter = 0;
    ll k = 0;
    ll k_saved = 0;
    ll hits = 0;
    ll sampled_ratio = 0;

    long double *norm_factor;

    ll *neo_previous;
    int *hist_WL;
    int *hist_E_selected;

    ll *neo_previous_root = 0;
    ll accept_counter_root = 0;
    ll reject_counter_root = 0;
    ll k_root = 0;
    ll k_saved_root = 0;
    ll hits_root = 0;
    ll sampled_ratio_root = 0;

    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // if ((rank == root && size == 1) || (rank == root && REP % (size - 1) != 0))
    // {
    //     cout << "It is required to have at least two processes running the algorithm. rank = 0 is the root and r != the workers" << endl;
    //     cout << "Or the number of worker processes might not be acceptable for this REP value (n_workers % REP != 0)." << endl;
    //     return 1;
    // }

    // Set a seed for the random number generator
    uint64_t seed;
    if (SEED == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }
    else
        seed = SEED;
    xorshift64s_state state = {.a = seed * rank};

    // NN file
    std::ifstream neighbour_tables_file(NEIGH_FILE);
    int *NN_table = new int[N_SPINS * NN];
    string line;

    if (neighbour_tables_file.is_open())
    {
        int i = 0;
        while (std::getline(neighbour_tables_file, line))
        {
            vector<string> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                NN_table[i++] = std::stold(a.at(idx));                
        }
        neighbour_tables_file.close();
    }
    else
    {
        cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << endl;
        return 1;
    }

    // Energies, magnetizations
    create_vector(energies, - max_E, max_E, 4);
    create_vector(magnetizations, - max_M, max_M, 2);

    ll REP = atoi(argv[2]);
    ll REP_per_worker = REP / (size - 1);

    if (rank == root)
    {
        int run = atoi(argv[1]);

        JDOS = new double[NE * NM];
        JDOS_frac = new double[NE * NM];

        for (int i = 0; i < NE * NM; i++)
        {
            JDOS[i] = 0;
            JDOS_frac[i] = 0;
        }

        JDOS[0] = 1;
        JDOS[NM - 1] = 1;
        JDOS_frac[0] = 1;
        JDOS_frac[NM - 1] = 1;

        // Normalization file
        std::ifstream norm_factor_file(NORM_FILE);
        string line;
        long double *norm_factor = new long double[NM];

        if (norm_factor_file.is_open()) 
        {
            for (int i = 0; std::getline(norm_factor_file, line); i++)
                norm_factor[i] = std::stold(line);
            norm_factor_file.close();
        }
        else 
        {
            cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << endl;
            return 1;
        }

        // Start measuring time
        vector<string> console_log;
        vector<string> data;
        time_t now;
        string t;
        
        now = time(0);
        t = ctime(&now); t.pop_back();

        string console_output = "run: " + std::to_string(run) + " | L: " + std::to_string(L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | workers: " + std::to_string(size - 1) + " (+1 root) | REP per worker " + std::to_string(REP_per_worker) + " | Lattice " + lattice;
        console_log.push_back(console_output);

        cout << endl;
        cout << console_output << endl;
        cout << "Starting time: " << t << endl << endl;

        auto method_start = std::chrono::steady_clock::now();

        for (q = 0; q < q_max - 1; q++)
        {
            MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(JDOS, NE * NM, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
            
            auto q_start = std::chrono::steady_clock::now();

            neo_previous_root = new ll[NE * NE];
            neo_previous = new ll[NE * NE];
            for (int i = 0; i < NE * NE; i++)
            {
                neo_previous_root[i] = 0;
                neo_previous[i] = 0;
            }
            
            accept_counter_root = 0;
            reject_counter_root = 0;
            k_root = 0;
            k_saved_root = 0;
            hits_root = 0;
            sampled_ratio_root = 0;

            MPI_Status statuses;
            double shuffle_time = 0;

            int repeat[size] = {0};
            int max_hits = 0;
            
            for (int i = 0; i < size - 1; i++)
            {
                MPI_Recv(neo_previous, NE * NE, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &statuses);
                for (int j = 0; j < NE * NE; j++)
                    neo_previous_root[j] += neo_previous[j];
               
                MPI_Recv(&accept_counter, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                accept_counter_root += accept_counter;

                MPI_Recv(&reject_counter, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                reject_counter_root += reject_counter;

                MPI_Recv(&k, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                k_root += k;

                MPI_Recv(&k_saved, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                k_saved_root += k_saved;

                MPI_Recv(&sampled_ratio, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sampled_ratio_root += sampled_ratio;

                MPI_Recv(&hits, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, statuses.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                hits_root += hits;

                if (hits > max_hits) max_hits = hits;
                if (hits < max_hits) cout << "<<< ERROR: One walker detected less E values >>>" << endl;

                if (statuses.MPI_SOURCE == root + 1)
                    MPI_Recv(&shuffle_time, 1, MPI_DOUBLE_PRECISION, root + 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            for (int i = 0; i < NE * NE; i++)
            {
                if (neo_previous_root[i] != 0)
                {
                    int idx_old = i / NE;
                    int idx_new = i % NE;

                    double sum_neo_previous_old = 0;
                    for (int j = 0; j < NE; j++)
                        sum_neo_previous_old += neo_previous_root[idx_old * NE + j];

                    JDOS_frac[idx_new * NM + (q + 1)] += JDOS_frac[idx_old * NM + q] * neo_previous_root[idx_old * NE + idx_new] / sum_neo_previous_old;
                }
            }

            for (int i = 0; i < NE; i++)
                JDOS[i * NM + (q + 1)] = JDOS_frac[i * NM + (q + 1)] * norm_factor[q + 1];

            auto q_end = std::chrono::steady_clock::now();
            double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

            now = time(0);
            t = ctime(&now); t.pop_back();

            string console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max - 2) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits_root / (size - 1)) + " | q_time/E: " + std::to_string((size - 1) * q_time / hits_root) + "s | shuffle time: " + std::to_string(shuffle_time) + "s";
            string data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits_root / (size - 1)) + " " + std::to_string((size - 1) * q_time / hits_root) +
            + " " + std::to_string(k_root) + " " + std::to_string(accept_counter_root) + " " + std::to_string(reject_counter_root) + " " + std::to_string(sampled_ratio_root / (size - 1));
            
            console_log.push_back(console_output);
            data.push_back(data_line);

            cout << console_output << endl;
        }

        // Stop mesuring time
        auto method_end = std::chrono::steady_clock::now();
        double runtime = (double) (std::chrono::duration_cast<std::chrono::microseconds> (method_end - method_start).count()) * pow(10, -6);
        now = time(0);
        t = ctime(&now); t.pop_back();

        std::cout << std::endl;
        std::cout << "Runtime: " << std::setw(8) << runtime << " seconds." << std::endl;
        std::cout << "Simulation ended at: " << t << std::endl;
        cout << SAVE_DIR(lattice, REP) << endl;
        // Write JDOS to file
        std::ofstream file1((string) SAVE_DIR(lattice, REP) + std::to_string(run) + "_" + SAVE_FILE(REP, skip) + ".txt");
        for (int i = 0; i < NE; i++) 
        {
            for (int j = 0; j < NM; j++) 
                file1 << JDOS[i * NM + j] << " ";
            file1 << "\n";
        }
        file1.close();

        std::ofstream file2((string) SAVE_DIR(lattice, REP) + std::to_string(run) + "_" + SAVE_FILE(REP, skip) + "_data.txt");
        file2 << "q q_max q_time hits q_time/hits k accept_counter reject_counter\n"; 
        for (int i = 0; i < data.size(); i++)
            file2 << data[i] << "\n";
        file2 << runtime << "\n";
        file2.close();

        std::ofstream file3((string) SAVE_DIR(lattice, REP) + std::to_string(run) + "_" + SAVE_FILE(REP, skip) + "_console_logs.txt");
        for (int i = 0; i < console_log.size(); i++)
            file3 << console_log.at(i) << "\n";
        file3.close();
        

        delete[] JDOS, JDOS_frac, norm_factor, neo_previous_root, neo_previous;
    }
    else
    {   do
        {
            MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
            if (q == 0)
            {
                JDOS = new double[NE * NM];
                neo_previous = new ll[NE * NE];
                hist_WL = new int[NE];
                hist_E_selected = new int[NE];

                neo_previous_root = new ll[NE * NE];
            }
            MPI_Bcast(JDOS, NE * NM, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
            
            for (int i = 0; i < NE * NE; i++)
            {
                neo_previous_root[i] = 0;
                neo_previous[i] = 0;
            }

            for (int i = 0; i < NE; i++)
            {
                hist_WL[i] = 0;
                hist_E_selected[i] = 0;
            }
            
            array<double, NE> WL_log_DOS;
            for (int i = 0; i < WL_log_DOS.size(); i++)
                WL_log_DOS[i] = log(JDOS[i * NM + q]);
            
            array<int, N_SPINS> spins_WL; spins_WL.fill(1);
            int E_WL_old = - max_E;

            vector<int> pos;
            vector<int> neg;
            for (int i = 0; i < N_SPINS; i++)
                pos.push_back(i);
            
            double shuffle_time = 0;

            if (q >= 1)
            {
                for (int idx = 0; idx < q; idx++)
                {
                    int idx_pos = xorshift64s(&state) % pos.size();
                    int flipped_pos = pos.at(idx_pos);

                    pos.erase(pos.begin() + idx_pos);
                    neg.push_back(flipped_pos);

                    spins_WL[flipped_pos] = - 1;

                    int sum_nei = 0;
                    for (int i = 0; i < NN; i++)
                        sum_nei += spins_WL[NN_table[flipped_pos * NN + i]];

                    int delta_E = - J * spins_WL[flipped_pos] * sum_nei;
                    E_WL_old += 2 * delta_E;
                }

                auto shuffle_start = std::chrono::high_resolution_clock::now();
                shuffle(spins_WL, REP, pos, neg, E_WL_old, NN_table, NN, J);
                auto shuffle_end = std::chrono::high_resolution_clock::now();

                shuffle_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (shuffle_end - shuffle_start).count()) * pow(10, -6);
            }

            int idx_E_WL_old = binary_search(energies, E_WL_old);
            int E_old = E_WL_old;
            int idx_E_old = idx_E_WL_old;
            
            fss_scan(J, NE, pos, spins_WL, neo_previous, NN_table, NN, E_old, idx_E_old, energies);
            
            hist_WL[idx_E_WL_old]++;
            hist_E_selected[idx_E_WL_old]++;
            
            accept_counter = 1;
            reject_counter = 0;

            k_saved = 1;
            k = 2;

            while (min_hist(hist_E_selected, NE) < REP_per_worker)
            {
                int idx_pos = xorshift64s(&state) % pos.size();
                int flipped_pos = pos.at(idx_pos);
                pos.erase(pos.begin() + idx_pos);
                neg.push_back(flipped_pos);

                spins_WL[flipped_pos] = - 1;

                int sum_nei = 0;
                for (int i = 0; i < NN; i++)
                    sum_nei += spins_WL[NN_table[flipped_pos * NN + i]];

                int delta_E = - J * spins_WL[flipped_pos] * sum_nei;
                int E_WL_new = E_WL_old + 2 * delta_E;
                
                int idx_neg = xorshift64s(&state) % neg.size();
                int flipped_neg = neg.at(idx_neg);
                neg.erase(neg.begin() + idx_neg);
                pos.push_back(flipped_neg);

                spins_WL[flipped_neg] = 1;
                
                sum_nei = 0;
                for (int i = 0; i < NN; i++)
                    sum_nei += spins_WL[NN_table[flipped_neg * NN + i]];

                delta_E = - J * spins_WL[flipped_neg] * sum_nei;
                E_WL_new += 2 * delta_E;
                int idx_E_WL_new = binary_search(energies, E_WL_new);

                double ratio = exp(WL_log_DOS[idx_E_WL_old] - WL_log_DOS[idx_E_WL_new]);

                if (ratio >= 1 || (ratio * 10000) > ((double) (xorshift64s(&state) % 10000)) || hist_WL[idx_E_WL_new] == 0)
                {
                    E_WL_old = E_WL_new;
                    idx_E_WL_old = idx_E_WL_new;
                    hist_WL[idx_E_WL_old]++;

                    accept_counter++;

                    if (k >= k_saved + skip && hist_E_selected[idx_E_WL_old] <= REP_per_worker || hist_WL[idx_E_WL_new] == 0) 
                    {
                        hist_E_selected[idx_E_WL_old]++;
                        k_saved = k;
                        
                        int E_old = E_WL_old;
                        int idx_E_old = idx_E_WL_old;

                        fss_scan(J, NE, pos, spins_WL, neo_previous, NN_table, NN, E_old, idx_E_old, energies);
                    }
                }
                else
                {
                    spins_WL[flipped_neg] = - 1;
                    spins_WL[flipped_pos] = 1;
                    
                    if (flipped_neg != flipped_pos)
                    {
                        pos.pop_back();
                        neg.pop_back();

                        pos.push_back(flipped_pos);
                        neg.push_back(flipped_neg);
                    }
                    
                    hist_WL[idx_E_WL_old]++;

                    reject_counter++;

                    if (k >= k_saved + skip && hist_E_selected[idx_E_WL_old] <= REP_per_worker) 
                    {
                        hist_E_selected[idx_E_WL_old]++;
                        k_saved = k;

                        int E_old = E_WL_old;
                        int idx_E_old = idx_E_WL_old;

                        fss_scan(J, NE, pos, spins_WL, neo_previous, NN_table, NN, E_old, idx_E_old, energies);
                    }
                }

                k++;
            }

            hits = 0;
            for (int i = 0; i < NE; i++)
                if (JDOS[i * NM + q] > 0)
                    hits++;

            sampled_ratio = REP_per_worker * hits * skip / k;

            MPI_Send(neo_previous, NE * NE, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&accept_counter, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&reject_counter, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&k, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&k_saved, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&sampled_ratio, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            MPI_Send(&hits, 1, MPI_LONG_LONG_INT, root, rank, MPI_COMM_WORLD);
            
            if (rank == root + 1)
                MPI_Send(&shuffle_time, 1, MPI_DOUBLE_PRECISION, root, 99, MPI_COMM_WORLD);
        } while (q != q_max - 2);
    }

    if (rank != root)
        delete[] JDOS, neo_previous, hist_WL, hist_E_selected, neo_previous_root;
    
    MPI_Finalize();

    return 0;
}


