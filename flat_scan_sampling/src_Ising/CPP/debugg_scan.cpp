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
#include <random>
#include <unordered_set>

using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;

#define SEED 0
#define ll long long

#define LATTICE_NUM 2
#define DIM_NUM 2
#define S 1/2
#define SZ 2 * (S + 1)

#define L 4
#define J 1

#define SAVE_DIR "./Data/"

#define N_SPINS L * L * L
#define NN 6
#define DIM "3D"
#define LATTICE "SC"

#define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#define NEIGH_FILE "./neighbour_tables/neighbour_table_3D_SC_6NN_L" + std::to_string(L) + ".txt"
#define SAVE_FILE(REP, skip) "JDOS_FSS_Ising_3D_SC_L" + std::to_string(L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip)

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
/* These real versions are due to Isaku Wada, 2002/01/09 added */

const int q_max = N_SPINS / 2 + 1;      // Max index of the magnetization computed (N_SPINS / 2 + 1 for half of the JDOS)
const int skip = N_SPINS;               // Scan a configuration each skip times to get more random results
const ll REP = pow(10, 4);              // Number of sampled configurations in one energy

const int max_E = (1.0 / 2.0) * NN * N_SPINS;
const int max_M = N_SPINS;

const int NE = 1 + (max_E / 2);         // Number of allowed energies
const int NM = N_SPINS + 1;             // Number of allowed magnetizations

std::unordered_set<int> pickSet(int n, int k, std::mt19937& gen)
{
    std::uniform_int_distribution<> dis(1, n);
    std::unordered_set<int> elems;

    while (elems.size() < k) {
        elems.insert(dis(gen));
    }

    return elems;
}

std::vector<int> pick(int n, int k) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::unordered_set<int> elems = pickSet(n, k, gen);

    // ok, now we have a set of k elements. but now
    // it's in a [unknown] deterministic order.
    // so we have to shuffle it:

    std::vector<int> result(elems.begin(), elems.end());
    std::shuffle(result.begin(), result.end(), gen);
    return result;
}



int main(int argc, char **argv)
{
    /*
RandStream.setGlobalStream( RandStream.create('mt19937ar','seed',100) );
pos = [1, 2, 3, 14, 11, 6, 10, 7];

for i = 1:length(pos)
    a(i) = randsample(pos,1);
end

for i = 1:100
    a(i) = pos(rem(fix(rand * 10000), length(pos)) + 1);
end
    */

    // genrand_res53();
    init_genrand(100);
    
    vector<int> pos = {1, 2, 3, 14, 11, 6, 10, 7};

    for (int i = 0; i < pos.size(); i++)
        cout << pick(10, 1).at(0) << endl;



    // for (int i = 0; i < 10; i++)
    //     cout << (int) (genrand_res53() * 10000) << endl;

    for (int i = 0; i < pos.size(); i++)
        cout << pos.at((int) (genrand_res53() * 10000) % pos.size()) << " ";
    cout << endl;

    return 0;
/*
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


    array<int, N_SPINS> spins;
    for (int i = 0; i < N_SPINS; i++)
    {
        spins[i] = 1;
        if (genrand_res53() < 0.75)
            spins[i] = -1;
    }

    vector<int> pos;
    for (int i = 0; i < N_SPINS; i++)
        if (spins[i] == 1) pos.push_back(i);

    int E_conf = 0; 
    for (int i = 0; i < N_SPINS; i++)
    {
        int sum_nei = 0;
        for (int a = 0; a < NN; a++)
            sum_nei += spins[NN_table[i * NN + a]];
        
        E_conf += - J * sum_nei * spins[i];
    }
    E_conf /= 2;
    int idx_E = binary_search(energies, E_conf);

    int M_conf = 0;
    for (int i = 0; i < N_SPINS; i++)
        M_conf += spins[i];
    int q = binary_search(magnetizations, M_conf);

    for (int i = 0; i < NE * NE; i++)
        neo_previous[i] = 0;

    fss_scan(J, NE, pos, spins, neo_previous, NN_table, NN, E_conf, idx_E, energies);
    
    for (int i = 0; i < NE; i++)
    {
        for (int j = 0; j < NE; j++)
            cout << neo_previous[i * NE + j] << " ";
        cout << endl;
    }
    

    for (int i = 0; i < NE * NE; i++)
        if (neo_previous[i] != 0)
            cout << i << endl;
    
    delete[] JDOS, JDOS_frac, norm_factor, NN_table;
    delete[] neo_previous, hist_E_selected, hist_WL;

    return 0;
    */
}

