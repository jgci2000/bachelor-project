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
#include "fss_functions.h"

using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;

#define ll long long

//////////////////////////////////////////////////////////////////////////////////////
//
//  Flat Scan Sampling for the 2D Ising spin S Model - v10
//  João Inácio, Jan. 21, 2021
//
//

// Seed for the rng. If SEED == 0, then the seed is random. 
#define SEED 0

// LATTICE -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex 
#define LATTICE 1
// DIM -> 1 - 2D; 2 - 3D 
#define DIM 1
// Spin-S particles 
#define S 1
// Number of spin projections 
#define SZ ((2 * S) + 1)

// Lattice size 
#define L 3
// Ineteraction strength 
#define J 1



// Saving directory for the JDOS file 
#define SAVE_DIR "./Data/"

#if DIM == 1 && LATTICE == 1
    // Number of particles 
#   define N_SPINS L * L
    // Number of nearest neughbours 
#   define NN 4
#   define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#   define NEIGH_FILE "./neighbour_tables/neighbour_table_2D_SS_4NN_L" + std::to_string(L) + ".txt"
#   define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_2D_SS_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#endif

#if DIM == 2 && LATTICE > 1 && LATTICE  < 7
#   if LATTICE == 2
        // Number of particles 
#       define N_SPINS L * L * L
        // Number of nearest neughbours 
#       define NN 6
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_SC_6NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_SC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE == 3
        // Number of particles 
#       define N_SPINS 2 * L * L * L
        // Number of nearest neughbours 
#       define NN 8
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_BCC_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_BCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE == 4
        // Number of particles 
#       define N_SPINS 4 * L * L * L
        // Number of nearest neughbours 
#       define NN 12
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_FCC_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_FCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE == 5
        // Number of particles 
#       define N_SPINS 2 * L * L * L
        // Number of nearest neughbours 
#       define NN 12
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_HCP_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_HCP_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE == 6
        // Number of particles 
#       define N_SPINS L * L * L
        // Number of nearest neughbours 
#       define NN 8
#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_Hex_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_Hex_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif

#endif

// Sum of Npos file
#define SUM_NPOS_FILE "./sum_npos/sum_configs_Npos" + std::to_string(SZ) + "_N_atm" + std::to_string(N_SPINS) + ".txt"

const int N = 2 * S * N_SPINS + 1;
// Max index of the magnetization computed (N_SPINS / 2 + 1 for half of the JDOS)
#if N % 2 == 0
const int q_max = N / 2 - 1;      
#else
const int q_max = (N + 1) / 2 - 1;
#endif

const int skip = N_SPINS;               // Scan a configuration each skip times to get more random results
const ll REP = pow(10, 3);              // Number of sampled configurations in one energy

const int max_E = 4 * S * S * N_SPINS * NN / 2;
const int max_M = 2 * S * N_SPINS;

const int NE = 1 + (max_E / 2);         // Number of allowed energies
const int NM = max_M + 1;               // Number of allowed magnetizations

// Random number generator
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
    // Seed for the RNG
    uint64_t seed;
    if (SEED == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }
    else
        seed = SEED;
    xorshift64s_state state = {.a = seed};
    
    // Normalization, NN file and sum  Npos
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

    std::ifstream neighbour_tables_file(NEIGH_FILE);
    int *NN_table = new int[N_SPINS * NN];

    if (neighbour_tables_file.is_open())
    {
        int i = 0;
        while (std::getline(neighbour_tables_file, line))
        {
            vector<int> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                NN_table[i++] = a.at(idx);
        }
        neighbour_tables_file.close();
    }
    else
    {
        cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << endl;
        return 1;
    }

    std::ifstream sum_npos_file(SUM_NPOS_FILE);
    int **sum_Npos;
    vector<int> line_size_sum_Npos;
    int size = 0;

    if (sum_npos_file.is_open())
    {
        int i = 0;

        std::getline(sum_npos_file, line);
        line_size_sum_Npos = split(line, ' ');

        sum_Npos = new int*[line_size_sum_Npos.size()];

        while (std::getline(sum_npos_file, line))
        {
            sum_Npos[i] = new int[line_size_sum_Npos.at(i) * SZ];
            vector<int> a = split(line, ' ');

            for (int k = 0; k < line_size_sum_Npos.at(i) * SZ; k++)
                sum_Npos[i][k] = a.at(k);
            i++;
        }
        sum_npos_file.close();
    }
    else
    {
        cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << endl;
        return 1;
    }
    
    // Energies, magnetizatios and JDOS
    array<int, NE> energies;
    array<int, NM> magnetizations;
    create_vector(energies, - max_E, max_E, 4);
    create_vector(magnetizations, - max_M, max_M, 2);

    array<int, SZ> Z_spin;
    create_vector(Z_spin, - 2 * S, 2 * S, 2);

    long double **JDOS_advance = new long double*[NM];
    
    for (int i = 0; i < NM; i++)
        JDOS_advance[i] = new long double[line_size_sum_Npos.at(i) * NE];
    
    for (int i = 0; i < NM; i++)
        for (int j = 0; j < NE * line_size_sum_Npos.at(i); j++)
            JDOS_advance[i][j] = 0;
    JDOS_advance[0][0] = 1;
    JDOS_advance[NM - 1][0] = 1;
    
    long double *JDOS = new long double[NE * NM];
    
    for (int i = 0; i < NE * NM; i++)
        JDOS[i] = 0;
    JDOS[0] = 1;
    JDOS[NM - 1] = 1;

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;
    
    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "L: " + std::to_string(L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip);
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();
    
    // New RPS algorithm
    
    /* Scan norm q1 */
    int *SPM = new int[N_SPINS];
    for (int i = 0; i < N_SPINS; i++) 
        SPM[i] = 0;
    
    int prev_E_old_idx = 0;
    int prev_Npos_sum_conf_old_idx = 0;

    for (int x = 0; x < SZ - 1; x++)
    {
        vector<int> flip_list;
        for (int i = 0; i < N_SPINS; i++)
            if (SPM[i] <= SZ - x - 2)
                flip_list.push_back(i);
        
        for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
        {
            int SPM_tmp[N_SPINS];
            for (int i = 0; i < N_SPINS; i++)
                SPM_tmp[i] = SPM[i];

            int E_tmp1 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp1 += - pow(Z_spin[SPM[flip_list.at(flip_idx)]], 2);
            
            int E_tmp2 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp2 += - Z_spin[SPM[flip_list.at(flip_idx)] + x + 1] * Z_spin[SPM[flip_list.at(flip_idx)]]; 
            
            int E_tmp3 =  E_tmp2 - E_tmp1 - max_E;
            SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x + 1;

            int counter[SZ] = {0};
            for (int i = 0; i < N_SPINS; i++)
                for (int j = 0; j < SZ; j++)
                    if (SPM_tmp[i] == j)
                        counter[j]++;
            
            int counter2 = 0;
            int idx_sum_Npos_vec;
            for (idx_sum_Npos_vec = 0; idx_sum_Npos_vec < line_size_sum_Npos.at(x + 1); idx_sum_Npos_vec += SZ)
            {
                for (int i = 0; i < SZ; i++)
                    if (counter[i] == sum_Npos[x + 1][i + idx_sum_Npos_vec])
                        counter2++;

                if (counter2 == SZ)
                    break;
                else
                    counter2 = 0;
            }
            idx_sum_Npos_vec /= SZ;

            int idx_E_tmp3 = binary_search(energies, E_tmp3);
            JDOS_advance[x + 1][idx_sum_Npos_vec * NE + idx_E_tmp3] += JDOS_advance[0][prev_Npos_sum_conf_old_idx * NE + prev_E_old_idx];
        }
    }

    double sum_JDOS_advance[NE] = {0};
    double sum_sum_JDOS_advance = 0;
    for (int i = 0; i < NE; i++)
    {
        for (int j = 0; j < line_size_sum_Npos.at(1); j++)
            sum_JDOS_advance[i] += JDOS_advance[1][j * NE + i];
        sum_sum_JDOS_advance += sum_JDOS_advance[i];
    }

    for (int i = 0; i < NE; i++)
        JDOS[i * NM + 1] = sum_JDOS_advance[i] * norm_factor[1] / sum_sum_JDOS_advance;

    // for (int idx = 0; idx < NM; idx++)
    // {
    //     for (int i = 0; i < NE; i++)
    //         {
    //             for (int j = 0; j < line_size_sum_Npos.at(idx); j++)
    //                 cout << JDOS_advance[idx][j * NE + i] << " ";
    //             cout << endl;
    //         }
    //     cout << endl;
    // }

    // for (int i = 0; i < N_SPINS; i++)
    //     cout << SPM[i] << " ";
    // cout << endl;
    
    // return 0;

    now = time(0);
    t = ctime(&now); t.pop_back();
    console_output = t + " | q: " + std::to_string(0) + "/" + std::to_string(q_max);
    console_log.push_back(console_output);
    
    cout << console_output << endl;

    ll *hist;
    ll *hist_E_selected;

    for (int q = 1; q <= q_max; q++) 
    {
        auto q_start = std::chrono::steady_clock::now();
        
        hist = new ll[NE * line_size_sum_Npos.at(q)];
        hist_E_selected = new ll[NE * line_size_sum_Npos.at(q)];

        for (int i = 0; i < NE * line_size_sum_Npos.at(q); i++)
        {
            hist[i] = 0;
            hist_E_selected[i] = 0;
        }

        /* random_spin_config_at_q */
        array<int, N_SPINS> spins_vector; spins_vector.fill(Z_spin[0]);
        for (int i = 0; i < N_SPINS; i++)
            SPM[i] = 0;
        
        int E_config = - max_E;

        for (int idx = 1; idx <= q; idx++)
        {
            vector<int> flip_list;
            for (int i = 0; i < N_SPINS; i++)
                if (SPM[i] < SZ - 1)
                    flip_list.push_back(i);

            // HERE
            int flip_list_idx = flip_list.at(xorshift64s(&state) % flip_list.size());
            // do flip_list_idx = flip_list.at(xorshift64s(&state) % flip_list.size());
            // while (SPM[flip_list_idx] >= SZ);

            // if (q == 6)
            // {
            //     for (int i = 0; i < N_SPINS; i++)
            //         cout << SPM[i] << " ";
            //     cout << endl;

            //     for (int i = 0; i < flip_list.size(); i++)
            //         cout << flip_list.at(i) << " ";
            //     cout << endl;

            //     cout << flip_list_idx << endl;
            // }
            
            

            int E_z_old_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_z_old_tmp += - spins_vector[flip_list_idx] * spins_vector[NN_table[flip_list_idx * NN + a]];
            
            SPM[flip_list_idx]++;
            spins_vector[flip_list_idx] = Z_spin.at(SPM[flip_list_idx]);

            int E_z_new_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_z_new_tmp += - spins_vector[flip_list_idx] * spins_vector[NN_table[flip_list_idx * NN + a]];

            E_config += E_z_new_tmp - E_z_old_tmp;
        }

        // if (q == 6)
        //     return 0;

        // for (int i = 0; i < N_SPINS; i++)
        //     cout << spins_vector[i] << " ";
        // cout << endl;

        // for (int i = 0; i < N_SPINS; i++)
        //     cout << SPM[i] << " ";
        // cout << endl; 

        int counter[SZ] = {0};
        for (int i = 0; i < N_SPINS; i++)
            for (int j = 0; j < SZ; j++)
                if (SPM[i] == j)
                    counter[j]++;

        int counter2 = 0;
        int idx_sum_Npos_vec;
        for (idx_sum_Npos_vec = 0; idx_sum_Npos_vec < line_size_sum_Npos.at(q); idx_sum_Npos_vec += SZ)
        {
            for (int i = 0; i < SZ; i++)
                if (counter[i] == sum_Npos[q][i + idx_sum_Npos_vec])
                    counter2++;

            if (counter2 == SZ)
                break;
            else
                counter2 = 0;
        }
        idx_sum_Npos_vec /= SZ;

        int idx_E_config = binary_search(energies, E_config);
        hist[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec]++;
        hist_E_selected[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec]++;

        // cout << E_config << endl;
        // for (int i = 0; i < NE; i++)
        // {
        //     for (int j = 0; j < line_size_sum_Npos.at(q); j++)
        //         cout << hist[i * line_size_sum_Npos.at(q) + j] << " ";
        //     cout << endl;
        // }

        /* scan_norm_correct */
        int prev_idx_E_config = idx_E_config;
        int prev_Npos_sum_conf_old_idx = idx_sum_Npos_vec;

        for (int x = 0; x < SZ - 1; x++)
        {
            vector<int> flip_list;
            for (int i = 0; i < N_SPINS; i++)
                if (SPM[i] < SZ - x - 1)
                    flip_list.push_back(i);
            
            // for (int i = 0; i < flip_list.size(); i++)
            //     cout << flip_list.at(i) << " ";
            // cout << endl;

            for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
            {
                int SPM_tmp[N_SPINS];
                for (int i = 0; i < N_SPINS; i++)
                    SPM_tmp[i] = SPM[i];

                int E_tmp1 = 0;
                for (int a = 0; a < NN; a++)
                    E_tmp1 += - spins_vector[flip_list.at(flip_idx)] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];
                
                int E_tmp2 = 0;
                for (int a = 0; a < NN; a++)
                    E_tmp2 += - Z_spin[SPM[flip_list.at(flip_idx)] + x + 1] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];                  
                
                int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x + 1;

                int counter[SZ] = {0};
                for (int i = 0; i < N_SPINS; i++)
                    for (int j = 0; j < SZ; j++)
                        if (SPM_tmp[i] == j)
                            counter[j]++;

                int counter2 = 0;
                int idx_sum_Npos_vec;
                for (idx_sum_Npos_vec = 0; idx_sum_Npos_vec < line_size_sum_Npos.at(q + x + 1); idx_sum_Npos_vec += SZ)
                {
                    for (int i = 0; i < SZ; i++)
                        if (counter[i] == sum_Npos[q + x + 1][i + idx_sum_Npos_vec])
                            counter2++;

                    if (counter2 == SZ)
                        break;
                    else
                        counter2 = 0;
                }
                idx_sum_Npos_vec /= SZ;

                int idx_E_tmp3 = binary_search(energies, E_tmp3);
                JDOS_advance[q + x + 1][idx_sum_Npos_vec * NE + idx_E_tmp3] += JDOS_advance[q][prev_Npos_sum_conf_old_idx * NE + prev_idx_E_config] / REP;
            }
        }
        
        idx_sum_Npos_vec = prev_Npos_sum_conf_old_idx;
        idx_E_config = prev_idx_E_config;
        ll k = 1;

        // for (int idx = 0; idx <= 5; idx++)
        // {
        //     for (int i = 0; i < line_size_sum_Npos.at(idx); i++)
        //     {
        //         for (int j = 0; j < NE; j++)
        //             cout << JDOS_advance[idx][i * NE + j] << " ";
        //         cout << endl;
        //     }
        //     cout << endl;
        // }

        // return 0;

        int  c = 0;
        while (min_hist(hist_E_selected, NE * line_size_sum_Npos.at(q)) < REP)
        {
            // cout << "start ->" << endl;
            // cout << "spins: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << spins_vector[i] << " ";
            // cout << endl;

            // cout << "SPM: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << SPM[i] << " ";
            // cout << endl;

            // cout << "E: " << E_config << endl;
            // cout << "idx_sum_Npos: " << idx_sum_Npos_vec << endl;

            /* rw_step_at_q */
            array<int, N_SPINS> spins_vector_new = spins_vector;
            int E_config_old = E_config;
            int idx_E_config_old = idx_E_config;
            int idx_sum_Npos_vec_old = idx_sum_Npos_vec;

            int SPM_new[N_SPINS];
            for (int i = 0; i < N_SPINS; i++)
                SPM_new[i] = SPM[i];

            int flipped_idx_1 = xorshift64s(&state) % N_SPINS;
            int flipped_pos_start = SPM[flipped_idx_1];
            int E_old_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_old_tmp += - spins_vector_new[flipped_idx_1] * spins_vector_new[NN_table[flipped_idx_1 * NN + a]];
            
            vector<int> SPM_end_list;
            for (int i = 0; i < SZ; i++)
                if (i != flipped_pos_start)
                    SPM_end_list.push_back(i);
            
            int SPM_end = SPM_end_list.at(xorshift64s(&state) % SPM_end_list.size());

            int SPM_start = SPM[flipped_idx_1];
            // int SPM_end;
            // do SPM_end = xorshift64s(&state) % SZ;
            // while (SPM_end == flipped_idx_1);

            spins_vector_new[flipped_idx_1] = Z_spin[SPM_end];
            SPM_new[flipped_idx_1] = SPM_end;

            int E_new_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_new_tmp += - spins_vector_new[flipped_idx_1] * spins_vector_new[NN_table[flipped_idx_1 * NN + a]];
            
            int SPM_diff = SPM_end - SPM_start;
            E_config = E_config_old - E_old_tmp + E_new_tmp;

            // cout << "after 1st flip ->" << endl;
            // cout << "spins: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << spins_vector_new[i] << " ";
            // cout << endl;

            // cout << "SPM: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << SPM_new[i] << " ";
            // cout << endl;

            // cout << "E: " << E_config << endl;
            // cout << "idx_sum_Npos: " << idx_sum_Npos_vec << endl;
            // cout << "diff: " << SPM_diff << endl;

            vector<int> flipped_list;
            for (int i = 0; i < N_SPINS; i++)
                if (SPM_new[i] - SPM_diff >= 0 && SPM_new[i] - SPM_diff < SZ)
                    flipped_list.push_back(i);
            
            int flipped_idx_2 = flipped_list.at(0);
            if (flipped_list.size() > 1)
                flipped_idx_2 = flipped_list.at(xorshift64s(&state) % flipped_list.size());
            // cout << flipped_idx_2 << endl;

            // int flipped_idx_2;
            // do flipped_idx_2 = xorshift64s(&state) % SZ;
            // while (SPM_new[flipped_idx_2] - SPM_diff < 0 || SPM_new[flipped_idx_2] - SPM_diff > SZ);

            E_old_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_old_tmp += - spins_vector_new[flipped_idx_2] * spins_vector_new[NN_table[flipped_idx_2 * NN + a]];

            SPM_new[flipped_idx_2] += - SPM_diff;
            // cout << SPM_new[flipped_idx_2] << endl;
            spins_vector_new[flipped_idx_2] = Z_spin[SPM_new[flipped_idx_2]];

            E_new_tmp = 0;
            for (int a = 0; a < NN; a++)
                E_new_tmp += - spins_vector_new[flipped_idx_2] * spins_vector_new[NN_table[flipped_idx_2 * NN + a]];

            E_config += - E_old_tmp + E_new_tmp;
            idx_E_config = binary_search(energies, E_config);

            int counter[SZ] = {0};
            for (int i = 0; i < N_SPINS; i++)
                for (int j = 0; j < SZ; j++)
                    if (SPM_new[i] == j)
                        counter[j]++;

            int counter2 = 0;
            for (idx_sum_Npos_vec = 0; idx_sum_Npos_vec < line_size_sum_Npos.at(q); idx_sum_Npos_vec += SZ)
            {
                for (int i = 0; i < SZ; i++)
                    if (counter[i] == sum_Npos[q][i + idx_sum_Npos_vec])
                        counter2++;

                if (counter2 == SZ)
                    break;
                else
                    counter2 = 0;
            }
            idx_sum_Npos_vec /= SZ;

            // if (q == 4)
            // {
            //     cout << "-----------------------------------------------------------" << endl;

            //     cout << "BEFORE: " << endl;
            //     for (int i = 0; i < N_SPINS; i++)
            //         cout << spins_vector[i] << " ";
            //     cout << endl;

            //     for (int i = 0; i < N_SPINS; i++)
            //         cout << SPM[i] << " ";
            //     cout << endl;

            //     cout << "E: " << E_config_old << "; idx: " << idx_E_config_old << endl;
            //     cout << "idx_sum_Npos: " << idx_sum_Npos_vec_old << endl;


            //     cout << "AFTER: " << endl;
            //     for (int i = 0; i < N_SPINS; i++)
            //         cout << spins_vector_new[i] << " ";
            //     cout << endl;

            //     for (int i = 0; i < N_SPINS; i++)
            //         cout << SPM_new[i] << " ";
            //     cout << endl;

            //     cout << "E: " << E_config << "; idx: " << idx_E_config << endl;
            //     cout << "idx_sum_Npos: " << idx_sum_Npos_vec << endl;

            //     for (int i = 0; i < NE; i++)
            //         cout << energies[i] << " ";
            //     cout << endl;
                
            //     cout << "-----------------------------------------------------------" << endl;
            //     c++;
            //     if (c == 100000)
            //     return 0;
            // }

            // cout << "end- >" << endl;
            // cout << "spins: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << spins_vector_new[i] << " ";
            // cout << endl;

            // cout << "SPM: " << endl;
            // for (int i = 0; i < N_SPINS; i++)
            //     cout << SPM_new[i] << " ";
            // cout << endl;

            // cout << "E: " << E_config << endl;
            // cout << "idx_sum_Npos: " << idx_sum_Npos_vec << endl;
            
            /* WL_rw_criteria */
            double ratio = JDOS_advance[q][idx_sum_Npos_vec * NE + idx_E_config] / JDOS_advance[q][idx_sum_Npos_vec_old * NE + idx_E_config_old];
            if (ratio >= 1 || xorshift64s(&state) % 10000 < ratio * 10000 || hist_E_selected[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec] == 0)
            {
                spins_vector = spins_vector_new;
                for (int i = 0; i < N_SPINS; i++)
                    SPM[i] = SPM_new[i];
            }
            else
            {
                E_config = E_config_old;
                idx_sum_Npos_vec = idx_sum_Npos_vec_old;

                idx_E_config = idx_E_config_old;
                // idx_sum_Npos_vec = idx_sum_Npos_vec_old;
            }

            hist[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec]++;

            if (hist_E_selected[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec] < REP && k % skip == 0)
            {
                /* scan_norm_correct */
                prev_idx_E_config = idx_E_config;
                prev_Npos_sum_conf_old_idx = idx_sum_Npos_vec;

                for (int x = 0; x < SZ - 1; x++)
                {
                    vector<int> flip_list;
                    for (int i = 0; i < N_SPINS; i++)
                        if (SPM[i] <= SZ - x - 2)
                            flip_list.push_back(i);

                    for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
                    {
                        int SPM_tmp[N_SPINS];
                        for (int i = 0; i < N_SPINS; i++)
                            SPM_tmp[i] = SPM[i];

                        int E_tmp1 = 0;
                        for (int a = 0; a < NN; a++)
                            E_tmp1 += - spins_vector[flip_list.at(flip_idx)] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];
                        
                        int E_tmp2 = 0;
                        for (int a = 0; a < NN; a++)
                            E_tmp2 += - Z_spin[SPM[flip_list.at(flip_idx)] + x + 1] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];                  
                        
                        int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                        SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x + 1;

                        int counter[SZ] = {0};
                        for (int i = 0; i < N_SPINS; i++)
                            for (int j = 0; j < SZ; j++)
                                if (SPM_tmp[i] == j)
                                    counter[j]++;

                        int counter2 = 0;
                        int idx_sum_Npos_vec;
                        for (idx_sum_Npos_vec = 0; idx_sum_Npos_vec < line_size_sum_Npos.at(q + x + 1); idx_sum_Npos_vec += SZ)
                        {
                            for (int i = 0; i < SZ; i++)
                                if (counter[i] == sum_Npos[q + x + 1][i + idx_sum_Npos_vec])
                                    counter2++;

                            if (counter2 == SZ)
                                break;
                            else
                                counter2 = 0;
                        }
                        idx_sum_Npos_vec /= SZ;

                        int idx_E_tmp3 = binary_search(energies, E_tmp3);
                        JDOS_advance[q + x + 1][idx_sum_Npos_vec * NE + idx_E_tmp3] += JDOS_advance[q][prev_Npos_sum_conf_old_idx * NE + prev_idx_E_config] / REP;
                    }
                }

                idx_sum_Npos_vec = prev_Npos_sum_conf_old_idx;
                idx_E_config = prev_idx_E_config;

                hist_E_selected[idx_E_config * line_size_sum_Npos.at(q) + idx_sum_Npos_vec]++;
            }

            k++;
        }

        
        

        // return 0;

        double sum_JDOS_advance[NE] = {0};
        double sum_sum_JDOS_advance = 0;
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < line_size_sum_Npos.at(q + 1); j++)
                sum_JDOS_advance[i] += JDOS_advance[q + 1][j * NE + i];
            sum_sum_JDOS_advance += sum_JDOS_advance[i];
        }

        for (int j = 0; j < line_size_sum_Npos.at(q + 1); j++)
            for (int i = 0; i < NE; i++)
                JDOS_advance[q + 1][j * NE + i] = JDOS_advance[q + 1][j * NE + i] * norm_factor[q + 1] / sum_sum_JDOS_advance;
        
        
        for (int i = 0; i < NE; i++)
            sum_JDOS_advance[i] = 0;
        // double sum_JDOS_advance[NE] = {0};
        sum_sum_JDOS_advance = 0;
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < line_size_sum_Npos.at(q + 1); j++)
                sum_JDOS_advance[i] += JDOS_advance[q + 1][j * NE + i];
            sum_sum_JDOS_advance += sum_JDOS_advance[i];
        }

        // if (q == 1)
        // {
        //     cout << sum_sum_JDOS_advance << endl;

        //     for (int i = 0; i < NE; i++)
        //         cout << sum_JDOS_advance[i] << " ";
        //     cout << endl;
        // } 

        for (int i = 0; i < NE; i++)
            JDOS[i * NM + q + 1] = sum_JDOS_advance[i] * norm_factor[q + 1] / sum_sum_JDOS_advance;

        int hits = 0;
        for (int j = 0; j < line_size_sum_Npos.at(q); j++)
            for (int i = 0; i < NE; i++)
                if (JDOS_advance[q][j * NE + i] > 0)
                    hits++;
        
        // for (int idx = 0; idx <= q_max + 1; idx++)
        // {
        //     for (int i = 0; i < line_size_sum_Npos.at(idx); i++)
        //     {
        //         for (int j = 0; j < NE; j++)
        //             cout << JDOS_advance[idx][i * NE + j] << " ";
        //         cout << endl;
        //     }
        //     cout << endl;
        // }

        // if (q == 6)
        // {
        //     for (int idx = 0; idx <= 5; idx++)
        //     {
        //         for (int i = 0; i < line_size_sum_Npos.at(idx); i++)
        //         {
        //             for (int j = 0; j < NE; j++)
        //                 cout << JDOS_advance[idx][i * NE + j] << " ";
        //             cout << endl;
        //         }
        //         cout << endl;
        //     }
        // }

        // if (q == 6)
        // {
        //     for (int i = 0; i < NE; i++)
        //     {
        //         for (int j = 0; j < NM; j++)
        //             cout << JDOS[i * NM + j] << " ";
        //         cout << endl;
        //     }

        //     cout << endl;

        //     for (int i = 0; i < NE; i++)
        //         cout << JDOS[i * NM + q + 1] << " ";
        //     cout << endl;

        //     return 0;
        // }
        

        // for (int i = 0; i < line_size_sum_Npos.at(2); i++)
        // {
        //     for (int j = 0; j < NE; j++)
        //         cout << JDOS_advance[2][i * NE + j] << " ";
        //     cout << endl;
        // }

        // cout << endl; 

        // for (int i = 0; i < NE; i++)
        // {
        //     for (int j = 0; j < line_size_sum_Npos.at(2); j++)
        //         cout << JDOS_advance[2][i * line_size_sum_Npos.at(2) + j] << " ";
        //     cout << endl;
        // }

        // for (int i = 0; i < N_SPINS; i++)
        //     cout << SPM[i] << " ";
        // cout << endl;

        // for (int i = 0; i < N_SPINS; i++)
        //     cout << spins_vector[i] << " "; 
        // cout << endl;

        // cout << "idx_sum_Npos_vec: " << idx_sum_Npos_vec << endl;

        delete[] hist, hist_E_selected;

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits); // + "s | shuffle time: " + std::to_string(shuffle_time) + "s";
        // data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        // + " " + std::to_string(k) + " " + std::to_string(accept_counter) + " " + std::to_string(reject_counter);
        
        console_log.push_back(console_output);
        // data.push_back(data_line);

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
    std::ofstream file1((string) SAVE_DIR + SAVE_FILE(REP, skip) + ".txt");
    for (int i = 0; i < NE; i++) 
    {
        for (int j = 0; j < NM; j++) 
            file1 << JDOS[i * NM + j] << " ";
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2((string) SAVE_DIR + SAVE_FILE(REP, skip) + "_data.txt");
    file2 << "q q_max q_time hits q_time/hits k accept_counter reject_counter\n"; 
    for (int i = 0; i < data.size(); i++)
        file2 << data[i] << "\n";
    file2 << runtime << "\n";
    file2.close();

    std::ofstream file3((string) SAVE_DIR + SAVE_FILE(REP, skip) + "_console_logs.txt");
    for (int i = 0; i < console_log.size(); i++)
        file3 << console_log.at(i) << "\n";
    file3.close();

    for (int i = 0; i < line_size_sum_Npos.size(); i++)
        delete[] sum_Npos[i];
    for (int i = 0; i < NM; i++)
        delete[] JDOS_advance[i];
    
    delete[] sum_Npos, JDOS_advance;
    delete[] JDOS, norm_factor, NN_table;
    delete[] SPM;
    
    return 0;
}


