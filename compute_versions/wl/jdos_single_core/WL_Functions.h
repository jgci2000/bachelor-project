//
// Header file for FSS method functions.
// João Inácio, Mar. 25th, 2021
//


#ifndef WL_FUNCTIONS_H
#define WL_FUNCTIONS_H

long long min_hist(long long *, int);

long double average_hist(long long *, int);

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

std::vector<std::string> split(const std::string& , char);



#endif // FSS_FUNCTIONS_H
