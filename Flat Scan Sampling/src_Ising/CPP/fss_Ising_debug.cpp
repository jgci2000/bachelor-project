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
//  Flat Scan Sampling for the 2D Ising Model - v10
//  João Inácio, Jan. 21, 2021
//
//

/* Seed for the rng. If SEED == 0, then the seed is random. */
#define SEED 0

/* LATTICE -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex */
#define LATTICE_NUM 1
/* DIM -> 1 - 2D; 2 - 3D */
#define DIM_NUM 1
/* S-Spins particles */
#define S 1/2
/* Number of spin projections */
#define SZ 2 * (S + 1)

/* Lattice size */
#define L 8
/* Ineteraction strength */
#define J 1



/* Saving directory for the JDOS file */
#define SAVE_DIR "./Data/"

#if DIM_NUM == 1 && LATTICE_NUM == 1
    /* Number of particles */
#   define N_SPINS L * L
    /* Number os nearest neughbours */
#   define NN 4
#   define DIM "2D"
#   define LATTICE "SS"

#   define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#   define NEIGH_FILE "./neighbour_tables/neighbour_table_2D_SS_4NN_L" + std::to_string(L) + ".txt"
#   define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_2D_SS_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#endif

#if DIM_NUM == 2 && LATTICE_NUM > 1 && LATTICE_NUM  < 7
#   if LATTICE_NUM == 2
        /* Number of particles */
#       define N_SPINS L * L * L
        /* Number os nearest neughbours */
#       define NN 6
#       define DIM "3D"
#       define LATTICE "SC"

#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_SC_6NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_SC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE_NUM == 3
        /* Number of particles */
#       define N_SPINS 2 * L * L * L
        /* Number os nearest neughbours */
#       define NN 8
#       define DIM "3D"
#       define LATTICE "BCC"

#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_BCC_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_BCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE_NUM == 4
        /* Number of particles */
#       define N_SPINS 4 * L * L * L
        /* Number os nearest neughbours */
#       define NN 12
#       define DIM "3D"
#       define LATTICE "FCC"

#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_FCC_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_FCC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE_NUM == 5
        /* Number of particles */
#       define N_SPINS 2 * L * L * L
        /* Number os nearest neughbours */
#       define NN 12
#       define DIM "3D"
#       define LATTICE "HCP"

#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_HCP_12NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_HCP_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif
#   if LATTICE_NUM == 6
        /* Number of particles */
#       define N_SPINS L * L * L
        /* Number os nearest neughbours */
#       define NN 8
#       define DIM "3D"
#       define LATTICE "Hex"

#       define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#       define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_Hex_8NN_L" + std::to_string(L) + ".txt"
#       define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_Hex_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)
#   endif

#endif

const int q_max = N_SPINS / 2 + 1;      // Max index of the magnetization computed (N_SPINS / 2 + 1 for half of the JDOS)
const int skip = N_SPINS;               // Scan a configuration each skip times to get more random results
const ll REP = pow(10, 3);              // Number of sampled configurations in one energy

const int max_E = (1.0 / 2.0) * NN * N_SPINS;
const int max_M = N_SPINS;

const int NE = 1 + (max_E / 2);         // Number of allowed energies
const int NM = N_SPINS + 1;             // Number of allowed magnetizations

/* Random number generator */
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

#include <stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 


int main(int argc, char **argv) 
{
    init_genrand(100);
    /* Seed for the RNG */
    uint64_t seed;
    if (SEED == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }
    else
        seed = SEED;
    xorshift64s_state state = {.a = seed};
    
    // Normalization and NN file
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

    // Energies, magnetizatios and JDOS
    array<int, NE> energies;
    array<int, NM> magnetizations;
    create_vector(energies, - max_E, max_E, 4);
    create_vector(magnetizations, - max_M, max_M, 2);
    
    long double *JDOS = new long double[NE * NM];
    long double *JDOS_frac = new long double[NE * NM];

    for (int i = 0; i < NE * NM; i++)
    {
        JDOS[i] = 0;
        JDOS_frac[i] = 0;
    }
    JDOS[0] = 1;
    JDOS[NM - 1] = 1;
    JDOS_frac[0] = 1;
    JDOS_frac[NM - 1] = 1;

    ll *neo_previous = new ll[NE * NE];
    ll *hist_WL = new ll[NE];
    ll *hist_E_selected = new ll[NE];

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;
    
    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "L: " + std::to_string(L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | dim: " + DIM + " | lattie: " + LATTICE;
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();
    
    // New RPS algorithm
    for (int q = 0; q < q_max - 1; q++) 
    {
        int count = 0;
        auto q_start = std::chrono::steady_clock::now();
        
        for (int i = 0; i < NE * NE; i++)
            neo_previous[i] = 0;

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
                int idx_pos = (int) (10000 * genrand_res53()) % pos.size();
                count++;
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

            // auto shuffle_start = std::chrono::high_resolution_clock::now();
            // shuffle(spins_WL, REP, pos, neg, E_WL_old, NN_table, NN, J);
            // auto shuffle_end = std::chrono::high_resolution_clock::now();

            // shuffle_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (shuffle_end - shuffle_start).count()) * pow(10, -6);
        }

        int idx_E_WL_old = binary_search(energies, E_WL_old);
        int E_old = E_WL_old;
        int idx_E_old = idx_E_WL_old;

        fss_scan(J, NE, pos, spins_WL, neo_previous, NN_table, NN, E_old, idx_E_old, energies);

        hist_WL[idx_E_WL_old]++;
        hist_E_selected[idx_E_WL_old]++;
        
        ll accept_counter = 1;
        ll reject_counter = 0;

        ll k_saved = 1;
        ll k = 2;
        
        while (min_hist(hist_E_selected, NE) < REP)
        {
            // if (q >= 5)
            //     cout << "min_H: " << min_hist(hist_E_selected, NE) << endl;
            int idx_pos = (int) (10000 * genrand_res53()) % pos.size();
            count++;
            int flipped_pos = pos.at(idx_pos);
            pos.erase(pos.begin() + idx_pos);
            neg.push_back(flipped_pos);

            spins_WL[flipped_pos] = - 1;

            int sum_nei = 0;
            for (int i = 0; i < NN; i++)
                sum_nei += spins_WL[NN_table[flipped_pos * NN + i]];

            int delta_E = - J * spins_WL[flipped_pos] * sum_nei;
            int E_WL_new = E_WL_old + 2 * delta_E;

            int idx_neg = (int) (10000 * genrand_res53()) % neg.size();
            count++;
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
            count++;
            if (ratio >= 1 || ratio > genrand_res53() || hist_E_selected[idx_E_WL_new] == 0)
            {
                E_WL_old = E_WL_new;
                idx_E_WL_old = idx_E_WL_new;
                hist_WL[idx_E_WL_old]++;

                accept_counter++;

                if (k >= k_saved + skip && hist_E_selected[idx_E_WL_old] < REP || hist_E_selected[idx_E_WL_new] == 0)
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

                if (k >= k_saved + skip && hist_E_selected[idx_E_WL_old] < REP)
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

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        string console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max - 2) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits) + "s | shuffle time: " + std::to_string(shuffle_time) + "s";
        string data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        + " " + std::to_string(k) + " " + std::to_string(accept_counter) + " " + std::to_string(reject_counter);
        
        console_log.push_back(console_output);
        data.push_back(data_line);

        cout << console_output << endl;
        cout << count << endl;
        cout << endl;

        // for (int i = 0; i < NE; i++)
        //     cout << hist_E_selected[i] << " ";
        // cout << endl;

        for (int i = 0; i < NE; i++)
            if (hist_WL[i] > 0)
                cout << hist_WL[i] << " ";
        cout << endl;

        // if (q == 15)
        // {
        //     for (int i = 0; i < NE / 4; i++)
        //         cout << JDOS[i * NM + q + 1] << endl;
        //     return 0;
        // }
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
    
    delete[] JDOS, JDOS_frac, norm_factor, NN_table;
    delete[] neo_previous, hist_E_selected, hist_WL;

    return 0;
}


