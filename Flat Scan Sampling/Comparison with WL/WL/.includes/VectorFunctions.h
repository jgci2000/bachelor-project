#ifndef _VECTOR_FUNCIIONS_H
#define _VECTOR_FUNCTIONS_H

#include <array>
#include <vector>

#include "VectorFunctions.cpp"

template<typename T, size_t N>
void create_vector(std::array<T, N> &, T, T, int);

template<typename T, size_t N>
void shift(std::array<T, N> &, std::array<T, N> &, int);

template<typename T, size_t N>
std::array<T, N> circshift(std::array<T, N> &in, int xshift, int yshift);

template<typename T, size_t N>
void print_1D(std::array<T, N> &);

template<typename T, size_t N>
void print_2D_square_vec(std::array<T, N> &);

template<typename T, size_t N>
std::vector<T> find(std::array<T, N> &, char, T);

template<typename T>
std::vector<T> find(T *, int, char, T);

template<typename T, size_t N>
int binary_search(std::array<T, N> &, T);

template<typename T>
int binary_search(T *arr, int size, T num);

#endif
