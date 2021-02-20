#ifndef _ISING_FUNCTIONS_H
#define _ISING_FUNCTIONS_H

#include "IsingFunctions.cpp"

template<typename T, size_t N>
void rand_spins2D(std::array<T, N> &);

template<typename T, size_t N1, size_t N2>
void comp_JDOS(std::array<T, N1> &, std::array<T, N2> &);

template<typename T, size_t N1, size_t N2>
int comp_energy2D(std::array<T, N1> &, std::array<T, N2> &, std::array<T, N2> &, const int);

template<typename T, size_t N>
int comp_magnetization2D(std::array<T, N> &);

template<typename T, size_t N1, size_t N2>
int comp_energy3D(std::array<T, N1> &, std::array<T, N2> &, std::array<T, N2> &, const int);

template<typename T, size_t N>
int comp_magnetization3D(std::array<T, N> &);

template<typename T, size_t N>
double average_hist(std::array<T, N> &);

template<typename T, size_t N>
int min_hist(std::array<T, N> &);

template<typename T>
T min_hist(T *hist, int size);

#endif