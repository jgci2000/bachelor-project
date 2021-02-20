#include <iostream>
#include <cmath>
#include <array>
#include <vector>

template<typename T, size_t N>
void print_1D(std::array<T, N> &array) 
{
    // Prints a 1D array.

    for (unsigned int i = 0; i < array.size(); i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
}

template<typename T, size_t N> 
void print_2D_square_vec(std::array<T, N> &array) 
{
    // Prints a 2D square vectorized array.

    int L = (int) sqrt(array.size());

    for (int i = 0; i < L; i++) 
    {
        for (int j = 0; j < L; j++)
            std::cout << array[i * L + j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename T, size_t N> 
void create_vector(std::array<T, N> &in, T init, T final, int step) 
{
    // Creates a vector "in" with values from "init" to "final" with a certain step.
    
    for (unsigned int i = 0; i < in.size(); i++)
        in[i] = init + i * step;
}

template<typename T, size_t N> 
void shift(std::array<T, N> &out, std::array<T, N> &in, int shft_amt) 
{
    // Shifts data in "in" by "shft_amt" and and records it in "out".

    for (unsigned int i = 0; i < in.size(); i++) 
    {
        int ii = (i + shft_amt) % in.size();
        if (ii < 0) ii = ii + in.size();
        out[ii] = in[i];
    }
}

template<typename T, size_t N>
std::array<T, N> circshift(std::array<T, N> &in, int xshift, int yshift)
{
    // Shift array circularly.

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

template<typename T, size_t N>
std::vector<T> find(std::array<T, N> &in, char op = '!', T num = 0)
{
    // Returns an array with the indices of nonzero items in the input.
    // If necessary, just add more options for the operation.

    std::vector<T> out;

    switch (op)
    {
    case '!':
        for (int i = 0; i < in.size(); i++) 
            if (in[i] != num)
                out.push_back(i);                
        break;
        
    case '=':
        for (int i = 0; i < in.size(); i++) 
            if (in[i] == num)
                out.push_back(i);
        break;
    }

    return out;
}

template<typename T>
std::vector<T> find(T *in, int size, char op = '!', T num = 0)
{
    // Returns an array with the indices of nonzero items in the input.
    // If necessary, just add more options for the operation.

    std::vector<T> out;

    switch (op)
    {
    case '!':
        for (int i = 0; i < size; i++) 
            if (in[i] != num)
                out.push_back(i);
        break;
        
    case '=':
        for (int i = 0; i < size; i++) 
            if (in[i] == num)
                out.push_back(i);
        break;
    }

    return out;
}

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


