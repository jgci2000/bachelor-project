#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "./.includes/VectorFunctions.h"
#include "./.includes/IsingFunctions.h"


using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;

//////////////////////////////////////////////////////////////////////////////////////
//
//  Wang-Landau Sampling for the 2D Ising Model
//  João Inácio, Dec. 26, 2020
//
//


#define WRITE_FILE  1                   // Enables writing the JDOS to a text file
#define SEED        0
#define lli         long long int

const int J = 1;                        // Interaction strength
const int L = 4;                        // Lattice size
const int N_SPINS = L * L;              // Number of particles
const int NN = 4;                       // Number of neighbours of each site

double f = exp(1);                      // Modification factor
// const double f_final = 1 + pow(10, -8);
const double flatness = 0.90;           // Flatness

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
    int exp_f_final = atoi(argv[1]);
    int run = atoi(argv[2]);

    double f_final = 1 + pow(10, - exp_f_final);

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
    
    // Fisrt configuration
    std::array<int, N_SPINS> spins_WL;
    for (int i = 0; i < spins_WL.size(); i++) 
    {
        if ((xorshift64s(&state) % 2) + 1 == 1) spins_WL[i] = + 1;
        else spins_WL[i] = - 1;
    }

    int E = 0, M = 0;

    for (int i = 0; i < N_SPINS; i++)
    {
        E += - (1.0 / 2.0) * J * spins_WL[i] * (spins_WL[nnxpos[i]] + spins_WL[nnxneg[i]] + spins_WL[nnypos[i]] + spins_WL[nnyneg[i]]);
        M += spins_WL[i];
    }

    int idx_E = binary_search(energies, E);
    int idx_M = binary_search(magnetizations, M);

    // JDOS and Histogram
    double *ln_JDOS = new double[NE * NM];
    double *hist = new double[NE * NM];

    for (int i = 0; i < NE * NM; i++)
    {
        ln_JDOS[i] = 0;
        hist[i] = 0;
    }

    int mc_sweeps = 0;

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;
    
    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "Run: " + std::to_string(run) + " | L: " + std::to_string(L) + " | f_final: " + std::to_string((int) - log10(f_final - 1)) + " | flatness: " + std::to_string((int) (flatness * 100));
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();
    std::chrono::_V2::steady_clock::time_point loop_start;

    while (f > f_final) 
    {
        if (mc_sweeps == 0)
            loop_start = std::chrono::steady_clock::now();

        for (int ni = 0; ni < N_SPINS; ni++) 
        {
            int i = (xorshift64s(&state) % N_SPINS);

            int new_E = E + 2 * J * spins_WL[i] * (spins_WL[nnxpos[i]] + spins_WL[nnxneg[i]] + spins_WL[nnypos[i]] + spins_WL[nnyneg[i]]);
            int new_M  = M - 2 * spins_WL[i];

            int idx_new_E = binary_search(energies, new_E);
            int idx_new_M = binary_search(magnetizations, new_M);
            
            double ratio = exp(ln_JDOS[idx_E * NM + idx_M] - ln_JDOS[idx_new_E * NM + idx_new_M]);

            if (ratio >= 1 || (ratio * 10000) > ((double) (xorshift64s(&state) % 10000)))
            {
                spins_WL[i] = - spins_WL[i];
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
            double avg_h = average_hist(hist, NE * NM);
            int min_h = min_hist(hist, NE * NM);

            if (min_h > avg_h * flatness) 
            {
                auto loop_end = std::chrono::steady_clock::now();
                double loop_dur = (double) (std::chrono::duration_cast<std::chrono::microseconds> (loop_end - loop_start).count()) * pow(10, -6);

                now = time(0);
                t = ctime(&now); t.pop_back();
                
                string console_output = t + " | f: " + std::to_string( - log10(f - 1)) + "/" + std::to_string((int) - log10(f_final - 1)) + " | Time elapsed: " + std::to_string(loop_dur) + "s";
                string data_line = std::to_string( - log10(f - 1)) + " " + std::to_string((int) - log10(f_final - 1)) + " " + std::to_string(loop_dur) + " " + std::to_string(mc_sweeps) + " " + std::to_string(min_h)
                + " " + std::to_string(avg_h);

                console_log.push_back(console_output);
                data.push_back(data_line);

                cout << console_output << endl;

                f = sqrt(f);
                mc_sweeps = 0;
                for (int i = 0; i < NE * NM; i++)
                    hist[i] = 0;
            }
        }
    }

    // Get JDOS
    double *JDOS = new double[NE * NM];
    
    for (int i = 0; i < NE * NM; i++)
        if (ln_JDOS[i] != 0) 
            JDOS[i] = exp(ln_JDOS[i] - ln_JDOS[0] + log(2)) / 2;    

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
        std::ofstream file1("./CPP Data/" + std::to_string(L) + "/" + std::to_string(exp_f_final) + "/" + std::to_string(run) + "_JDOS_2D_" + std::to_string(L) + "L_f" + std::to_string((int) - log10(f_final - 1)) + "_flatness" + std::to_string((int) (flatness * 100)) + "_WL_CPP.txt");
        for (int i = 0; i < NE; i++) 
        {
            for (int j = 0; j < NM; j++) 
                file1 << JDOS[i * NM + j] << " ";
            file1 << "\n";
        }
        file1.close();

        std::ofstream file2("./CPP Data/" + std::to_string(L) + "/" + std::to_string(exp_f_final) + "/" + std::to_string(run) + "_JDOS_2D_" + std::to_string(L) + "L_f" + std::to_string((int) - log10(f_final - 1)) + "_flatness" + std::to_string((int) (flatness * 100)) + "_WL_CPP_data.txt");
        file2 << "f f_max loop_dur mc_sweeps min_h avg_h\n"; 
        for (int i = 0; i < data.size(); i++)
            file2 << data[i] << "\n";
        file2 << runtime << "\n";
        file2.close();

        std::ofstream file3("./CPP Data/" + std::to_string(L) + "/" + std::to_string(exp_f_final) + "/" + std::to_string(run) + "_JDOS_2D_" + std::to_string(L) + "L_f" + std::to_string((int) - log10(f_final - 1)) + "_flatness" + std::to_string((int) (flatness * 100)) + "_WL_CPP_console_logs.txt");
        for (int i = 0; i < console_log.size(); i++)
            file3 << console_log.at(i) << "\n";
        file3.close();
    }

    delete[] JDOS, hist, ln_JDOS;

    return 0;
}


