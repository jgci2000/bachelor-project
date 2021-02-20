#include <iostream>
#include <cmath>
#include <array>
#include <climits>
#include <time.h>

template<typename T, size_t N>
void rand_spins(std::array<T, N> &spins) 
{
    // Creates a random +/- 1 vectorized 2D square spins matrix.

    srand((unsigned) time(NULL));

    for (unsigned int i = 0; i < spins.size(); i++) 
    {
        if ((rand() % 2) + 1 == 1) spins[i] = + 1;
        else spins[i] = - 1;
    }
}

template<typename T, size_t N>
void comp_JDOS(std::array<T, N> &lnJDOS, std::array<T, N> &JDOS) 
{
    // Computes the actual JDOS from an unormalized ln(JDOS) 2D vectorized matrix.

    double lnJDOS0 = lnJDOS[0];

    for (unsigned int i = 0; i < lnJDOS.size(); i++)
        if (lnJDOS[i] != 0) {
            lnJDOS[i] = lnJDOS[i] - lnJDOS0 + log(2);
            JDOS[i] = exp(lnJDOS[i]) / 2;
        }
}

template<typename T, size_t N1, size_t N2>
int comp_energy2D(std::array<T, N1> &spins, std::array<T, N2> &plus, std::array<T, N2> &minus, const int J) 
{
    // Computes the energy of a vectorized 2D Ising square lattice with periodic boudary conditions.

    int E = 0;
    int L = (int) sqrt(spins.size());

    for (int i = 0; i < L; i++) 
        for (int j = 0; j < L; j++) 
        {
            int s = spins[i * L + j];
            int sum_nei = spins[plus[i] * L + j] + spins[minus[i] * L + j] + spins[i * L + plus[j]] + spins[i * L + minus[j]];
            E += - J * s * sum_nei;
        }

    return E / 2;
}

template<typename T, size_t N>
int comp_magnetization2D(std::array<T, N> &spins) 
{
    // Computes the magnetization of a vectorized 2D Ising square lattice with periodic boudary conditions.

    int M = 0;

    for (unsigned int i = 0; i < spins.size(); i++)
        M += spins[i];
    
    return M;
}

template<typename T, size_t N1, size_t N2>
int comp_energy3D(std::array<T, N1> &spins, std::array<T, N2> &plus, std::array<T, N2> &minus, const int J) 
{
    // Computes the energy of a vectorized 3D Ising square lattice with periodic boudary conditions.

    int E = 0;
    int L = plus.size();
    int L2 = L * L;

    for (int i = 0; i < L; i++) 
        for (int j = 0; j < L; j++) 
            for (int k = 0; k < L; k++) 
            {
                int sum_nei = spins[j + plus[i] * L + k * L2] + spins[j + minus[i] * L + k * L2] + spins[plus[j] + i * L + k * L2] +
                 + spins[minus[j] + i * L + k * L2] + spins[j + i * L + plus[k] * L2] + spins[j + i * L + minus[k] * L2];
                E += - J * spins[j + i * L + k * L2] * sum_nei;
            }

    return E / 2;
}

template<typename T, size_t N>
int comp_magnetization3D(std::array<T, N> &spins) 
{
    // Computes the magnetization of a vectorized 3D Ising square lattice with periodic boudary conditions.

    int M = 0;

    for (unsigned int i = 0; i < spins.size(); i++)
        M += spins[i];
    
    return M;
}

template<typename T, size_t N>
double average_hist(std::array<T, N> &hist) 
{
    // Calculates the average of the histogram for the WL sampling.

    double sum = 0;
    int counter = 0;

    for (unsigned int i = 0; i < hist.size(); i++) 
        if (hist[i] != 0)
        {
            sum += hist[i];
            counter++;
        }
            
    return sum / counter;
}

template<typename T>
double average_hist(T *hist, int size) 
{
    // Calculates the average of the histogram for the WL sampling.

    double sum = 0;
    int counter = 0;

    for (unsigned int i = 0; i < size; i++) 
        if (hist[i] != 0)
        {
            sum += hist[i];
            counter++;
        }
            
    return sum / counter;
}

template<typename T, size_t N>
int min_hist(std::array<T, N> &hist) 
{
    // Finds the minimum value of the histogram for WL sampling.
    
    int min = INT_MAX;

    for (unsigned int i = 0; i < hist.size(); i++) 
        if (hist[i] != 0 && hist[i] < min)
            min = hist[i];

    return min;
}

template<typename T>
T min_hist(T *hist, int size) 
{
    // Finds the minimum value of the histogram for WL sampling.
    
    int min = INT_MAX;

    for (unsigned int i = 0; i < size; i++) 
        if (hist[i] != 0 && hist[i] < min)
            min = hist[i];

    return min;
}





