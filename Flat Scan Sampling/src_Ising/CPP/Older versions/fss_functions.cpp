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
void fss_scan(int J, int NE, std::vector<T> &pos, std::array<T, N1> &spins_WL, long long *neo_previous, int E_old, int idx_E_old, std::array<T, N2> &energies, std::array<T, N1> &nnxpos, std::array<T, N1> &nnxneg, std::array<T, N1> &nnypos, std::array<T, N1> &nnyneg)
{
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

        neo_previous[idx_E_old * NE + idx_E_new]++;
        spins_WL[flipped_pos_scan] = 1;
    }
}

