//
// Source file for FSS method functions.
// João Inácio, Mar. 25th, 2021
//


#include <climits>
#include <array>
#include <vector>
#include <string>

#include "WL_Functions.h"


// Finds the minimum value for the histogram, using C arrays.
long long min_hist(long long *hist, int size) 
{
    long long min = LONG_LONG_MAX;
    for (int i = 0; i < size; i++) 
        if (hist[i] != 0 && hist[i] < min)
            min = hist[i];
    return min;
}

// Finds the average of the histogram.
long double average_hist(long long *hist, int size)
{
    long double sum = 0;
    int nnz = 0;
    for (int i = 0; i < size; i++)
    {
        if (hist[i] != 0)
        {
            sum += hist[i];
            nnz++;
        }
    }
    return sum / nnz;
}

// Improved binary search algorithm, using C array.
int binary_search(int *arr, int size, int num)
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

// Improved binary search algorithm, using C++ std::vector.
int binary_search(std::vector<int> &arr, int num)
{
    int idx_low = -1; 
    int idx_high = arr.size(); 
    int idx_middle;
    while(idx_low + 1 != idx_high)
    {
        idx_middle = (idx_low + idx_high) / 2;
        if(arr.at(idx_middle) < num)
            idx_low = idx_middle;
        else
            idx_high = idx_middle;
    }
    return (idx_high >= arr.size() || arr.at(idx_high) != num) ? -1 : idx_high;
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
