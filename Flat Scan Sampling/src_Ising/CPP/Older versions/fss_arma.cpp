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
#include <armadillo>
#include "fss_functions.h"

using namespace arma;
using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::string;

//////////////////////////////////////////////////////////////////////////////////////
//
//  New Random Path Sampling for the 2D Ising Model - v9
//  João Inácio, Dec. 18, 2020
//
//

#define WRITE_FILE  1                   // Enables writing the JDOS to a text file
#define SEED        0
#define ll          long long

const int J = 1;                        // Interaction strength
const int L = 8;                        // Lattice size
const int N_SPINS = L * L;              // Number of particles
const int NN = 4;                       // Number of neighbours of each site

const int q_max = N_SPINS / 2 + 1;      // Max index of the magnetization computed (N_SPINS / 2 + 1 for half of the JDOS)
const int skip = N_SPINS;               // Scan a configuration each skip times to get more random results
const ll REP = pow(10, 3);   // Number of sampled configurations in one energy

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

template<class ty>
ty circshift(ty &in, int xshift, int yshift)
{
    ty out = in;
    for (int i = 0; i < in.n_rows; i++) 
    {
        int ii = (i + xshift) % in.n_rows;
        if (ii < 0) ii = in.n_rows + ii;
        for (int j = 0; j < in.n_cols; j++) 
        {
            int jj = (j + yshift) % in.n_cols;
            if (jj < 0) jj = in.n_cols + jj;
            out(ii, jj) = in(i, j);
        }
    }
    return out;
}

// Improved binary search algorithm, using Armadillo mat.
template<class ty, typename T>
int binary_search(ty &arr, T num)
{
    int idx_low = -1; 
    int idx_high = arr.n_elem; 
    int idx_middle;
    while(idx_low + 1 != idx_high)
    {
        idx_middle = (idx_low + idx_high) / 2;
        if(arr(idx_middle) < num)
            idx_low = idx_middle;
        else
            idx_high = idx_middle;
    }
    return (idx_high >= arr.size() || arr(idx_high) != num) ? -1 : idx_high;
}

template<class ty>
double min_hist(ty &hist)
{
    double min = INT_MAX;
    for (int i = 0; i < hist.n_elem; i++) 
        if (hist(i) != 0 && hist(i) < min)
            min = hist(i);
    return min;
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
    mat lattice(L, L);
    for (int i = 0; i < lattice.n_elem; i++)
        lattice(i) = i;
    
    mat nnxpos = circshift(lattice, 0, - 1);
    mat nnxneg = circshift(lattice, 0, 1);
    mat nnypos = circshift(lattice, 1, 0);
    mat nnyneg = circshift(lattice, - 1, 0);

    vec energies = arma::linspace(- max_E, max_E, NE);
    vec magnetizations = arma::linspace(- max_M, max_M, NM);
    
    mat JDOS(NE, NM, arma::fill::zeros);
    mat JDOS_frac(NE, NM, arma::fill::zeros);

    JDOS(0, 0) = 1;
    JDOS(0, JDOS.n_cols - 1) = 1;
    JDOS_frac(0) = 1;
    JDOS_frac(0, JDOS_frac.n_cols - 1) = 1;

    mat neo_previous(NE, NE, arma::fill::zeros);
    vec hist_WL(NE, arma::fill::zeros);
    vec hist_E_selected(NE, arma::fill::zeros);

    // Open the normalization factor
    std::ifstream norm_factor_file("./norm_factors/norm_factor_Ising_Natm_" + std::to_string(N_SPINS) + ".txt");
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
        return 1;
    }

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
    for (int q = 0; q < q_max - 1; q++) 
    {
        auto q_start = std::chrono::steady_clock::now();

        neo_previous.zeros();
        hist_WL.zeros();
        hist_E_selected.zeros();

        vec WL_log_DOS = arma::log(JDOS(arma::span::all, q));

        mat spins_WL(L, L, arma::fill::ones);
        int E_WL_old = - max_E;
        
        vector<int> pos;
        vector<int> neg;
        for (int i = 0; i < lattice.n_elem; i++)
            pos.push_back(lattice(i));

        if (q >= 1)
        {
            for (int idx = 0; idx < q; idx++)
            {
                int idx_pos = xorshift64s(&state) % pos.size();
                int flipped_pos = pos.at(idx_pos);

                pos.erase(pos.begin() + idx_pos);
                neg.push_back(flipped_pos);

                spins_WL(flipped_pos) = - 1;

                int delta_E = - J * spins_WL(flipped_pos) * 
                (spins_WL(nnxpos(flipped_pos)) + 
                spins_WL(nnxneg(flipped_pos)) + 
                spins_WL(nnypos(flipped_pos)) + 
                spins_WL(nnyneg(flipped_pos)));
                E_WL_old += 2 * delta_E;
            }
        }
        
        int idx_E_WL_old = binary_search(energies, E_WL_old);
        int E_old = E_WL_old;
        int idx_E_old = idx_E_WL_old;
        
        for (int i = 0; i < pos.size(); i++)
        {
            int flipped_pos_scan = pos.at(i);
            spins_WL(flipped_pos_scan) = - 1;

            int delta_E = - J * spins_WL(flipped_pos_scan) * 
                (spins_WL(nnxpos(flipped_pos_scan)) + 
                spins_WL(nnxneg(flipped_pos_scan)) + 
                spins_WL(nnypos(flipped_pos_scan)) + 
                spins_WL(nnyneg(flipped_pos_scan)));
            int E_new = E_old + 2 * delta_E;
            
            int idx_E_new = binary_search(energies, E_new);

            neo_previous(idx_E_old, idx_E_new)++;
            spins_WL(flipped_pos_scan) = 1;
        }
        
        hist_WL(idx_E_WL_old)++;
        hist_E_selected(idx_E_WL_old)++;
        
        ll accept_counter = 1;
        ll reject_counter = 0;

        ll k_saved = 1;
        ll k = 2;
        
        while (min_hist(hist_E_selected) < REP)
        {
            int idx_pos = xorshift64s(&state) % pos.size();
            int flipped_pos = pos.at(idx_pos);
            pos.erase(pos.begin() + idx_pos);
            neg.push_back(flipped_pos);

            spins_WL(flipped_pos) = - 1;

            int delta_E = - J * spins_WL(flipped_pos) * 
                (spins_WL(nnxpos(flipped_pos)) + 
                spins_WL(nnxneg(flipped_pos)) + 
                spins_WL(nnypos(flipped_pos)) + 
                spins_WL(nnyneg(flipped_pos)));
            int E_WL_new = E_WL_old + 2 * delta_E;
            
            int idx_neg = xorshift64s(&state) % neg.size();
            int flipped_neg = neg.at(idx_neg);
            neg.erase(neg.begin() + idx_neg);
            pos.push_back(flipped_neg);

            spins_WL(flipped_neg) = 1;
            
            delta_E = - J * spins_WL(flipped_neg) * 
                (spins_WL(nnxpos(flipped_neg)) + 
                spins_WL(nnxneg(flipped_neg)) + 
                spins_WL(nnypos(flipped_neg)) + 
                spins_WL(nnyneg(flipped_neg)));
            E_WL_new += 2 * delta_E;
            int idx_E_WL_new = binary_search(energies, E_WL_new);

            double ratio = exp(WL_log_DOS(idx_E_WL_old) - WL_log_DOS(idx_E_WL_new));

            if (ratio >= 1 || (ratio * 10000) > ((double) (xorshift64s(&state) % 10000)))
            {
                E_WL_old = E_WL_new;
                idx_E_WL_old = idx_E_WL_new;
                hist_WL(idx_E_WL_old)++;

                accept_counter++;

                if (k >= k_saved + skip && hist_E_selected(idx_E_WL_old) <= REP) 
                {
                    hist_E_selected(idx_E_WL_old)++;
                    k_saved = k;
                    
                    int E_old = E_WL_old;
                    int idx_E_old = idx_E_WL_old;

                    for (int i = 0; i < pos.size(); i++)
                    {
                        int flipped_pos_scan = pos.at(i);
                        spins_WL(flipped_pos_scan) = - 1;

                        int delta_E = - J * spins_WL(flipped_pos_scan) * 
                            (spins_WL(nnxpos(flipped_pos_scan)) + 
                            spins_WL(nnxneg(flipped_pos_scan)) + 
                            spins_WL(nnypos(flipped_pos_scan)) + 
                            spins_WL(nnyneg(flipped_pos_scan)));
                        int E_new = E_old + 2 * delta_E;

                        int idx_E_new = binary_search(energies, E_new);

                        neo_previous(idx_E_old, idx_E_new)++;
                        spins_WL(flipped_pos_scan) = 1;
                    }
                }
            }
            else
            {
                spins_WL(flipped_neg) = - 1;
                spins_WL(flipped_pos) = 1;
                
                if (flipped_neg != flipped_pos)
                {
                    pos.pop_back();
                    neg.pop_back();

                    pos.push_back(flipped_pos);
                    neg.push_back(flipped_neg);
                }
                
                hist_WL(idx_E_WL_old)++;

                reject_counter++;

                if (k >= k_saved + skip && hist_E_selected(idx_E_WL_old) <= REP) 
                {
                    hist_E_selected(idx_E_WL_old)++;
                    k_saved = k;

                    int E_old = E_WL_old;
                    int idx_E_old = idx_E_WL_old;

                    for (int i = 0; i < pos.size(); i++)
                    {
                        int flipped_pos_scan = pos.at(i);
                        spins_WL(flipped_pos_scan) = - 1;

                        int delta_E = - J * spins_WL(flipped_pos_scan) * 
                            (spins_WL(nnxpos(flipped_pos_scan)) + 
                            spins_WL(nnxneg(flipped_pos_scan)) + 
                            spins_WL(nnypos(flipped_pos_scan)) + 
                            spins_WL(nnyneg(flipped_pos_scan)));
                        int E_new = E_old + 2 * delta_E;

                        int idx_E_new = binary_search(energies, E_new);

                        neo_previous(idx_E_old, idx_E_new)++;
                        spins_WL(flipped_pos_scan) = 1;
                    }
                }
            }

            k++;
        }

        int hits = arma::sum(JDOS_frac(arma::span::all, q) > 0);

        for (int i = 0; i < NE * NE; i++)
        {
            if (neo_previous[i] != 0)
            {
                uvec idx_E = arma::ind2sub(size(neo_previous), i);
                JDOS_frac(idx_E(1), q + 1) += JDOS_frac(idx_E(0), q) * neo_previous(idx_E(0), idx_E(1)) / arma::sum(neo_previous(idx_E(0), arma::span::all));
            }
        }

        JDOS(arma::span::all, q + 1) = JDOS_frac(arma::span::all, q + 1) * norm_factor[q + 1];

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        string console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max - 2) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits) + "s";
        string data_line = std::to_string(q) + " " + std::to_string(q_max - 2) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        + " " + std::to_string(k) + " " + std::to_string(accept_counter) + " " + std::to_string(reject_counter);
        
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
        std::ofstream file1("./CPP Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP.txt");
        for (int i = 0; i < NE; i++) 
        {
            for (int j = 0; j < NM; j++) 
                file1 << JDOS(i, j) << " ";
            file1 << "\n";
        }
        file1.close();

        std::ofstream file2("./CPP Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP_data.txt");
        file2 << "q q_max q_time hits q_time/hits k accept_counter reject_counter\n"; 
        for (int i = 0; i < data.size(); i++)
            file2 << data.at(i) << "\n";
        file2 << runtime << "\n";
        file2.close();

        std::ofstream file3("./CPP Data/JDOS_2D_" + std::to_string(L) + "L_10E" + std::to_string((int) log10(REP)) + "_skip" + std::to_string(skip) + "_nRPS_CPP_console_logs.txt");
        for (int i = 0; i < console_log.size(); i++)
            file3 << console_log.at(i) << "\n";
        file3.close();
    }

    delete[] norm_factor;

    return 0;
}


