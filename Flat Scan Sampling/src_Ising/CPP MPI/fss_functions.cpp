#include <climits>
#include <array>
#include <vector>

#include "fss_functions.h"


// Finds the minimum value for the histogram, using C arrays.
template<typename T>
T min_hist(T *hist, int size) 
{
    int min = INT_MAX;
    for (int i = 0; i < size; i++) 
        if (hist[i] != 0 && hist[i] < min)
            min = hist[i];
    return min;
}

// Improved binary search algorithm, using C array.
template<typename T>
int binary_search(T *arr, int size, T num)
{
    int idx_low = -1; 
    int idx_high = size; 
    int idx_middle;
    while(idx_low + 1 != idx_high)
    {
        idx_middle = (idx_low + idx_high) / 2;
        if(arr[idx_middle] < num)
            idx_low = idx_middle;
        else
            idx_high = idx_middle;
    }
    return (idx_high >= size || arr[idx_high] != num) ? -1 : idx_high;
}

// Improved binary search algorithm, using C++ std::array.
template<typename T, size_t N>
int binary_search(std::array<T, N> &arr, T num)
{
    int idx_low = -1; 
    int idx_high = arr.size(); 
    int idx_middle;
    while(idx_low + 1 != idx_high)
    {
        idx_middle = (idx_low + idx_high) / 2;
        if(arr[idx_middle] < num)
            idx_low = idx_middle;
        else
            idx_high = idx_middle;
    }
    return (idx_high >= arr.size() || arr[idx_high] != num) ? -1 : idx_high;
}

// Creates a vector "in" with values from "init" to "final" with a certain step.
template<typename T, size_t N> 
void create_vector(std::array<T, N> &in, T init, T final, int step) 
{
    for (int i = 0; i < in.size(); i++)
        in[i] = init + i * step;
}

// Shift array circularly, using C++ std::array.
template<typename T, size_t N>
std::array<T, N> circshift(std::array<T, N> &in, int xshift, int yshift)
{
    std::array<T, N> out;
    int L = sqrt(in.size());
    for (int i =0; i < L; i++) 
    {
        int ii = (i + xshift) % L;
        if (ii < 0) ii = L + ii;
        for (int j = 0; j < L; j++) 
        {
            int jj = (j + yshift) % L;
            if (jj < 0) jj = L + jj;
            out[ii * L + jj] = in[i * L + j];
        }
    }
    return out;
}

// Scan for the Flat Scan Sampling.
template<typename T, size_t N1, size_t N2>
void fss_scan(int J, int NE, std::vector<T> &pos, std::array<T, N1> &spins_WL, long long *neo_previous, int *NN_table, int NN, int E_old, int idx_E_old, std::array<T, N2> &energies)
{
    for (int idx = 0; idx < pos.size(); idx++)
    {
        int flipped_pos_scan = pos.at(idx);
        spins_WL[flipped_pos_scan] = - 1;

        int sum_nei = 0;
        for (int i = 0; i < NN; i++)
            sum_nei += spins_WL[NN_table[flipped_pos_scan * NN + i]];

        int delta_E = - J * spins_WL[flipped_pos_scan] * sum_nei;
        int E_new = E_old + 2 * delta_E;

        int idx_E_new = binary_search(energies, E_new);

        neo_previous[idx_E_old * NE + idx_E_new]++;
        spins_WL[flipped_pos_scan] = 1;
    }
}

// Scan for the Flat Scan Sampling.
template<typename T, size_t N1, size_t N2>
void fss_scan_thread(int *ret, int J, int NE, std::vector<T> &pos, std::array<T, N1> &spins_WL, int *NN_table, int NN, int E_old, int idx_E_old, std::array<T, N2> &energies)
{
    for (int idx = 0; idx < pos.size(); idx++)
    {
        int flipped_pos_scan = pos.at(idx);
        spins_WL[flipped_pos_scan] = - 1;

        int sum_nei = 0;
        for (int i = 0; i < NN; i++)
            sum_nei += spins_WL[NN_table[flipped_pos_scan * NN + i]];

        int delta_E = - J * spins_WL[flipped_pos_scan] * sum_nei;
        int E_new = E_old + 2 * delta_E;

        int idx_E_new = binary_search(energies, E_new);

        ret[idx] = idx_E_old * NE + idx_E_new;
        spins_WL[flipped_pos_scan] = 1;
    }
}

// Shuffle spins.
template<typename T, size_t N>
void shuffle(std::array<T, N> &spins_WL, long long REP, std::vector<T> &pos, std::vector<T> &neg, int &E_WL_old, int *NN_table, int NN, int J)
{
    srand((unsigned) time(NULL));

    for (int i = 0; i < REP; i++)
    {
        int idx_pos = rand() % pos.size();
        int flipped_pos = pos.at(idx_pos);

        int idx_neg = rand() % neg.size();
        int flipped_neg = neg.at(idx_neg);

        if (flipped_neg == flipped_pos)
            continue;
        
        pos.erase(pos.begin() + idx_pos);
        neg.push_back(flipped_pos);

        spins_WL[flipped_pos] = - 1;

        int sum_nei = 0;
        for (int i = 0; i < NN; i++)
            sum_nei += spins_WL[NN_table[flipped_pos * NN + i]];

        int delta_E = - J * spins_WL[flipped_pos] * sum_nei;
        E_WL_old += 2 * delta_E;


        neg.erase(neg.begin() + idx_neg);
        pos.push_back(flipped_neg);

        spins_WL[flipped_neg] = 1;

        sum_nei = 0;
        for (int i = 0; i < NN; i++)
            sum_nei += spins_WL[NN_table[flipped_neg * NN + i]];

        delta_E = - J * spins_WL[flipped_neg] * sum_nei;
        E_WL_old += 2 * delta_E;
    }
}

// Split strings.
std::vector<std::string> split(const std::string& s, char seperator)
{
    std::vector<std::string> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(substring);

        prev_pos = ++pos;
    }

    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word

    return output;
}
