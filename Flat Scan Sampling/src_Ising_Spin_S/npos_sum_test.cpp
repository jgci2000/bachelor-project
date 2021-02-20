#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "fss_functions.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;

#define S 1
/* Number of spin projections */
#define SZ ((2 * S) + 1)

/* Lattice size */
#define L 4
/* Ineteraction strength */
#define J 1

#define N_SPINS L * L
#define NN 4

#   define NORM_FILE "./coefficients/coefficients_" + std::to_string(N_SPINS) + "d" + std::to_string(SZ) + ".txt"
#   define NEIGH_FILE "./neighbour_tables/neighbour_table_2D_SS_4NN_L" + std::to_string(L) + ".txt"


#define SUM_NPOS_FILE "./sum_npos/sum_configs_Npos" + std::to_string(SZ) + "_N_atm" + std::to_string(N_SPINS) + ".txt"



const int max_E = 4 * S * S * N_SPINS * NN / 2;
const int max_M = 2 * S * N_SPINS;

const int NE = 1 + (max_E / 2);         // Number of allowed energies
const int NM = max_M + 1;               // Number of allowed magnetizations

int main()
{
    cout << max_E << endl;
    cout << NE << endl;

    cout << max_M << endl;
    cout << NM << endl;

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

    // cout << "bytes allocated before sum_npos: " << sizeof(long double) * NM + sizeof(int) * N_SPINS * NN << endl;

    std::ifstream sum_npos_file(SUM_NPOS_FILE);
    int **sum_Npos;
    vector<int> line_size_sum_Npos;

    if (sum_npos_file.is_open())
    {
        int i = 0;

        std::getline(sum_npos_file, line);
        line_size_sum_Npos = split(line, ' ');
        for (int j = 0; j < line_size_sum_Npos.size(); j++)
            line_size_sum_Npos.at(j) *= SZ;
        
        sum_Npos = new int*[line_size_sum_Npos.size()];

        while (std::getline(sum_npos_file, line))
        {
            // cout << i << endl;
            // cout << "tryed allocating " << sizeof(int) * line_size_sum_Npos.at(i) << " bytes" << endl;
            sum_Npos[i] = new int[line_size_sum_Npos.at(i)];
            vector<int> a = split(line, ' ');

            for (int k = 0; k < line_size_sum_Npos.at(i); k++)
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

    // for (int i = 0; i < line_size_sum_Npos.size(); i++)
    //     cout << line_size_sum_Npos.at(i) << " ";
    // cout << endl;

    // for (int i = 0; i < line_size_sum_Npos.size(); i++)
    // {
    //     for (int j = 0; j < line_size_sum_Npos.at(i); j++)
    //         cout << sum_Npos[i][j] << " ";
    //     cout << endl;
    // }

    for (int i = 0; i < line_size_sum_Npos.size(); i++)
        delete[] sum_Npos[i];

    delete[] sum_Npos;
    delete[] norm_factor;
    delete[] NN_table;
    return 0;
}
