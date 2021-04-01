#include <iostream>
#include <cmath>
#include <array>

template<typename T, size_t N>
void print_spins(std::array<T, N> &spins) {
    // Prints a L*L pseudo matrix
    std::cout << "\tSpins Matrix:" <<std::endl;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++)
            std::cout << std::setw(3) << spins[i + j * L];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename T, size_t N>
void print_hist(std::array<T, N> &hist, std::array<T, N> &energies, std::array<T, N> &magnetizations) {
    std::cout << "\tHistogram:" << std::endl << "   ";

    for (int i = 0; i < cols; i++) {
        printf("%7d", magnetizations[i]);
    }
    std::cout << std::endl;

    for (int i = 0; i < rows; i++) {
        printf("%4d: ", energies[i]);
        for (int j = 0; j < cols; j++)
            printf("%7d", hist[i * cols + j]);
        std::cout << std::endl;
    }

    std::cout << std::endl;     
}

template<typename T, size_t N>
void print_JDOS(std::array<T, N> &JDOS, std::array<T, N> &energies, std::array<T, N> &magnetizations) {
    std::cout << "\tJDOS:" << std::endl << "   ";

    for (int i = 0; i < cols; i++) {
        printf("%8d", magnetizations[i]);
    }
    std::cout << std::endl;

    for (int i = 0; i < rows; i++) {
        printf("%4d: ", energies[i]);
        for (int j = 0; j < cols; j++)
            printf("%.2f\t", JDOS[i * cols + j]);
        std::cout << std::endl;
    }

    std::cout << std::endl;     
}

