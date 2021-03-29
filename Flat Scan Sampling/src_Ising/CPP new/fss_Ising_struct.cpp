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
#define SAVE_DIR    "./Data/"


struct Ising_struct
{
    int L;
    int N_atm;
    int NN;
    int J = 1;

    int NE;
    int NM;
    int max_E;
    int max_M;
    
    std::map<int, int> energies;
    std::map<int, int> magnetizations;

    int dim;
    std::string lattice;

    int *spins_vector;
    int *NN_table;
    long double *norm_factor;

    int E_config;
    int idx_E_config;

    int M_config;
    int idx_M_config;
};

void read_NN_talbe(Ising_struct *ising, std::string file_name)
{
    std::ifstream neighbour_tables_file(file_name);
    std::string line;

    if (neighbour_tables_file.is_open())
    {
        int i = 0;
        while (std::getline(neighbour_tables_file, line))
        {
            std::vector<std::string> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                ising->NN_table[i++] = std::stold(a.at(idx));
        }
        neighbour_tables_file.close();
    }
    else
        std::cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << std::endl;
}

void read_norm_factor(Ising_struct *ising, std::string file_name)
{
    std::ifstream norm_factor_file(file_name);
    std::string line;

    if (norm_factor_file.is_open()) 
    {
        for (int i = 0; std::getline(norm_factor_file, line); i++)
            ising->norm_factor[i] = std::stold(line);
        norm_factor_file.close();
    }
    else 
        std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
}

std::map<int, int> create_map(int init, int final, int step)
{
    std::map<int, int> out;
    int i = 0;
    while (init <= final)
    {
        out.insert(std::pair<int, int>(init, i));
        init += step;
        i++;
    }
    return out;
}

void set_E_config(Ising_struct *ising, int E_config)
{
    ising->E_config = E_config;
    ising->idx_E_config = ising->energies[E_config];
}

void get_system(Ising_struct *ising, int L, int lattice_num)
{
    switch (lattice_num)
    {
        case 1:
            ising->dim = 2;
            ising->lattice = "SS";
            ising->N_atm = L * L;
            ising->NN = 4;
            break;

        case 2:
            ising->dim = 3;
            ising->lattice = "SC";
            ising->N_atm = L * L * L;
            ising->NN = 6;
            break;

        case 3:
            ising->dim = 3;
            ising->lattice = "BCC";
            ising->N_atm = 2 * L * L * L; 
            ising->NN = 8;
            break;
        
        case 4:
            ising->dim = 3;
            ising->lattice = "FCC";
            ising->N_atm = 4 * L * L * L;
            ising->NN = 12;
            break;

        case 5: 
            ising->dim = 3;
            ising->lattice = "HCP";
            ising->N_atm = 2 * L * L * L;
            ising->NN = 12;
            break;
        case 6:
            ising->dim = 3;
            ising->lattice = "Hex";
            ising->N_atm = L * L * L;
            ising->NN = 8;
            break;

        default:
            std::cout << "Invalid lattice number." << std::endl;
            break;
    }
}


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
    Ising_struct *ising = new Ising_struct;
    ising->L = L_LATTICE;
    get_system(ising, L_LATTICE, LATTICE_NUM);

    ising->max_E = (1.0 / 2.0) * ising->NN * ising->N_atm;
    ising->max_M = ising->N_atm;

    ising->NE = 1 + (ising->max_E / 2);
    ising->NM = ising->N_atm + 1;

    ising->spins_vector = new int[ising->N_atm];
    ising->NN_table = new int[ising->N_atm * ising->NN];
    ising->norm_factor = new long double[ising->NM];

    ising->energies = create_map(- ising->max_E, ising->max_E, 4);
    ising->magnetizations = create_map(- ising->max_M, ising->max_M, 2);

    int q_max = (ising->NM + 1) / 2 - 2;

    if (ising->NM % 2 == 0)
        q_max = ising->NM / 2 - 3;

    int skip = ising->N_atm;
    ll REP = pow(10, 4);

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + std::to_string(ising->dim) + "D_" + ising->lattice + "_" + std::to_string(ising->NN) + "NN_L" + std::to_string(ising->L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + std::to_string(ising->N_atm) + "d2.txt";
    string save_file = "JDOS_FSS_Ising_" + std::to_string(ising->dim) + "D_" + ising->lattice + "_L" + std::to_string(ising->L) + "_REP_1E" + std::to_string((int) log10(REP)) + "_skip_" + std::to_string(skip);

    read_NN_talbe(ising, NN_table_file_name);
    read_norm_factor(ising, norm_factor_file_name);

    ld *JDOS = new ld[ising->NE * ising->NM];
    ll *hist = new ll[ising->NE];
    ll *hist_E_selected = new ll[ising->NE];

    for (int i = 0; i < ising->NE * ising->NM; i++)
        JDOS[i] = 0;
    JDOS[0] = 1;

    // Start measuring time
    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "L: " + std::to_string(ising->L) + " | REP: " + std::to_string(REP) + " | skip: " + std::to_string(skip) + " | dim: " + std::to_string(ising->dim) + "D | lattie: " + ising->lattice;
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // Flat Scan Sampling

    // Scan and Compute JDOS at q = 1    
    for (int flip_idx = 0; flip_idx < ising->N_atm; flip_idx++)
    {
        int E_tmp1 = 0;
        for (int a = 0; a < ising->NN; a++)
            E_tmp1 += - 1;
        
        int E_tmp2 = 0;
        for (int a = 0; a < ising->NN; a++)
            E_tmp2 += 1;
        
        int E_tmp3 = - ising->max_E - E_tmp1 + E_tmp2;
        int idx_E_tmp3 = ising->energies[E_tmp3];

        JDOS[idx_E_tmp3 * ising->NM + 1] += JDOS[0];
    }

    int sum_JDOS = 0;
    for (int i = 0; i < ising->NE; i++)
        if (JDOS[i * ising->NM + 1] > 0)
            sum_JDOS += JDOS[i * ising->NM + 1];

    for (int i = 0; i < ising->NE; i++)
        JDOS[i * ising->NM + 1] = JDOS[i * ising->NM + 1] * ising->norm_factor[1] / sum_JDOS;

    console_output = t + " | q: " + std::to_string(0) + "/" + std::to_string(q_max);
    console_log.push_back(console_output);

    cout << console_output << endl;

    // Main Loop    
    for (int q = 1; q <= q_max; q++)
    {
        auto q_start = std::chrono::steady_clock::now();

        for (int i = 0; i < ising->NE; i++)
        {
            hist[i] = 0; 
            hist_E_selected[i] = 0;
        }

        // Random config at q
        for (int i = 0; i < ising->N_atm; i++)
            ising->spins_vector[i] = 1;
        set_E_config(ising, - ising->max_E);

        array<vector<int>, 2> flip_list;
        for (int i = 0; i < ising->N_atm; i++)
            flip_list[0].push_back(i);

        int E_config = ising->E_config;
        for (int idx = 1; idx <= q; idx++)
        {
            int idx_tmp = rand_xoshiro256pp() % flip_list[0].size();
            int flipped_idx = flip_list[0].at(idx_tmp);
            ising->spins_vector[flipped_idx] = - 1;

            flip_list[1].push_back(flipped_idx);
            flip_list[0].erase(flip_list[0].begin() + idx_tmp);
            
            int delta_E = 0;
            for (int a = 0; a < ising->NN; a++)
                delta_E += - ising->spins_vector[flipped_idx] * ising->spins_vector[ising->NN_table[flipped_idx * ising->NN + a]];
            
            E_config += 2 * delta_E;
        }
        set_E_config(ising, E_config);

        // Update Histograms
        hist[ising->idx_E_config]++;
        hist_E_selected[ising->idx_E_config]++;

        // Scan the first config
        for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
        {
            int delta_E = 0;
            for (int a = 0; a < ising->NN; a++)
                delta_E += ising->spins_vector[flip_list[0].at(flip_idx)] * ising->spins_vector[ising->NN_table[flip_list[0].at(flip_idx) * ising->NN + a]];

            int E_tmp = ising->E_config + 2 * delta_E;
            int idx_E_tmp = ising->energies[E_tmp];

            JDOS[idx_E_tmp * ising->NM + q + 1] += JDOS[ising->idx_E_config * ising->NM + q] / REP;
        }

        ll k = 1;
        int *new_spins_vector = new int[ising->N_atm];

        // Where the magic happens
        while (min_hist(hist_E_selected, ising->NE) < REP)
        {
            // Get a new random condig at magnetization q
            for (int i = 0; i < ising->N_atm; i++)
                new_spins_vector[i] =  ising->spins_vector[i];
            int new_E_config = 0;
            int new_idx_E_config = 0;
            array<vector<int>, 2> new_flip_list = flip_list;

            // Flip a positive spin to a negative
            int idx_tmp1 = rand_xoshiro256pp() % new_flip_list[0].size();
            int flipped_idx1 = new_flip_list[0].at(idx_tmp1);
            new_spins_vector[flipped_idx1] = - 1;

            new_flip_list[1].push_back(flipped_idx1);
            new_flip_list[0].erase(new_flip_list[0].begin() + idx_tmp1);

            int delta_E = 0;
            for (int a = 0; a < ising->NN; a++)
                delta_E += - new_spins_vector[flipped_idx1] * new_spins_vector[ising->NN_table[flipped_idx1 * ising->NN + a]];
            new_E_config = ising->E_config + 2 * delta_E;

            // Flip a negative spin to a positive
            int idx_tmp2 = rand_xoshiro256pp() % new_flip_list[1].size();
            int flipped_idx2 = new_flip_list[1].at(idx_tmp2);
            new_spins_vector[flipped_idx2] = 1;

            new_flip_list[0].push_back(flipped_idx2);
            new_flip_list[1].erase(new_flip_list[1].begin() + idx_tmp2);

            delta_E = 0;
            for (int a = 0; a < ising->NN; a++)
                delta_E += - new_spins_vector[flipped_idx2] * new_spins_vector[ising->NN_table[flipped_idx2 * ising->NN + a]];
            new_E_config = new_E_config + 2 * delta_E;          

            // Wang Landau criteria
            new_idx_E_config = ising->energies[new_E_config];
            ld ratio = JDOS[ising->idx_E_config * ising->NM + q] / JDOS[new_idx_E_config * ising->NM + q];

            if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio)
            {
                for (int i = 0; i < ising->N_atm; i++)
                    ising->spins_vector[i] = new_spins_vector[i];

                set_E_config(ising, new_E_config);
                flip_list = new_flip_list;
            }

            hist[ising->idx_E_config]++;

            // Scan configuration
            if (hist_E_selected[ising->idx_E_config] < REP && k % skip == 0)
            {
                for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
                {
                    int delta_E = 0;
                    for (int a = 0; a < ising->NN; a++)
                        delta_E += ising->spins_vector[flip_list[0].at(flip_idx)] * ising->spins_vector[ising->NN_table[flip_list[0].at(flip_idx) * ising->NN + a]];

                    int E_tmp = ising->E_config + 2 * delta_E;
                    int idx_E_tmp = ising->energies[E_tmp];

                    JDOS[idx_E_tmp * ising->NM + q + 1] += JDOS[ising->idx_E_config * ising->NM + q] / REP;
                }

                hist_E_selected[ising->idx_E_config]++;
            }

            k++;
        }

        delete[] new_spins_vector;

        ld sum_JDOS = 0;
        for (int i = 0; i < ising->NE; i++)
            if (JDOS[i * ising->NM + q + 1] > 0)
                sum_JDOS += JDOS[i * ising->NM + q + 1];

                

        for (int i = 0; i < ising->NE; i++)
            JDOS[i * ising->NM + q + 1] = JDOS[i * ising->NM + q + 1] * ising->norm_factor[q + 1] / sum_JDOS;

        int hits = 0;
        for (int i = 0; i < ising->NE; i++)
            if (JDOS[i * ising->NM + q] > 0)
                hits++;

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        console_output = t + " | q: " + std::to_string(q) + "/" + std::to_string(q_max) + " | q_time: " + std::to_string(q_time) + "s | E: " + std::to_string(hits) + " | q_time/E: " + std::to_string(q_time / hits) + "s";
        string data_line = std::to_string(q) + " " + std::to_string(q_max) + " " + std::to_string(q_time) + " " + std::to_string(hits) + " " + std::to_string(q_time / hits) +
        + " " + std::to_string(k);
        
        console_log.push_back(console_output);
        data.push_back(data_line);

        cout << console_output << endl;

        // for (int i = 0; i < ising->NE; i++)
        //     if (JDOS[i * ising->NM + q + 1] > 0) cout << JDOS[i * ising->NM + q + 1] << " ";
        // cout << endl;

        // for (int i = 0; i < ising->NE; i++)
        //     if (hist[i] > 0) cout << hist[i] << " ";
        // cout << endl;
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
    for (int i = 0; i < ising->NE; i++) 
    {
        for (int j = 0; j < ising->NM; j++) 
            file1 << JDOS[i * ising->NM + j] << " ";
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2((string) SAVE_DIR + save_file + "_data.txt");
    file2 << "q q_max q_time hits q_time/hits k \n"; 
    for (int i = 0; i < data.size(); i++)
        file2 << data[i] << "\n";
    file2 << runtime << "\n";
    file2.close();

    std::ofstream file3((string) SAVE_DIR + save_file + "_console_logs.txt");
    for (int i = 0; i < console_log.size(); i++)
        file3 << console_log.at(i) << "\n";
    file3.close();
    
    delete[] JDOS, hist, hist_E_selected;
    delete[] ising->spins_vector, ising->NN_table, ising->norm_factor;

    return 0;
}
