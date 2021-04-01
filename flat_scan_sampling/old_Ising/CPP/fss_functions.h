#ifndef FSS_FUNCTIONS_H
#define FSS_FUNCTIONS_H

#include "fss_functions.cpp"

template<typename T>
T min_hist(T *, int);

template<typename T>
int binary_search(T *, int, T);

template<typename T, size_t N>
int binary_search(std::array<T, N> &, T);

template<typename T, size_t N>
void create_vector(std::array<T, N> &, T, T, int);

template<typename T, size_t N>
void shift(std::array<T, N> &, std::array<T, N> &, int);

template<typename T, size_t N>
std::array<T, N> circshift(std::array<T, N> &, int, int);

template<typename T, size_t N1, size_t N2>
void fss_scan(int, int, std::vector<T> &, std::array<T, N1> &, long long *, int *, int, int, int, std::array<T, N2> &);

template<typename T, size_t N>
void shuffle(std::array<T, N> &, long long, std::vector<T> &, std::vector<T> &, int &, int *, int, int);

std::vector<std::string> split(const std::string& , char);


#endif
