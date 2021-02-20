#ifndef _ISING_PRINTS_H
#define _ISING_PRINTS_H

#include "IsingPrints.cpp"

template<typename T, size_t N>
void print_spins(std::array<T, N> &);

template<typename T, size_t N>
void print_hist(std::array<T, N> &, std::array<T, N> &, std::array<T, N> &);

template<typename T, size_t N>
void print_JDOS(std::array<T, N> &, std::array<T, N> &, std::array<T, N> &);


#endif