//
// Header file for FSS method functions.
// João Inácio, Mar. 25th, 2021
//


#ifndef FSS_FUNCTIONS_H
#define FSS_FUNCTIONS_H

#include "Ising.h"

long long min_hist(long long *, int);

int binary_search(int *, int, int);

template<typename T, size_t N>
int binary_search(std::array<T, N> &, T);

int binary_search(std::vector<int> &, int);

template<typename T, size_t N>
void create_vector(std::array<T, N> &, T, T, int);

template<typename T, size_t N>
void shift(std::array<T, N> &, std::array<T, N> &, int);

template<typename T, size_t N>
std::array<T, N> circshift(std::array<T, N> &, int, int);

template<typename T, size_t N1, size_t N2>
void fss_scan(int, int, std::vector<T> &, std::array<T, N1> &, long long *, int *, int, int, int, std::array<T, N2> &);

void shuffle(Ising &, std::array<std::vector<int>, 2> &, int, long long, long double *);

std::vector<std::string> split(const std::string& , char);

#endif // FSS_FUNCTIONS_H
