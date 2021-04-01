#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <algorithm>

void print1D(double* arr,  int size) {
    for (int i = 0; i < size; i++)
        std::cout << arr[i] << " ";
    std::cout << std::endl;
}

void printSpins(int* spins, int L) {
    // Prints a L*L vectorized matrix
    std::cout << "\tSpins Matrix:" <<std::endl;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++)
            std::cout << std::setw(3) << spins[i + j * L];
        std::cout << std::endl;
    }
}

double* createVector(double init, double final, double step) {
    // Creates a vector like init:step:final in MATLAB
    int N = ((final - init) / step) + 1;
    double* out = new double[N];
    for (int i = 0; i < N; i++) 
        out[i] = init + i * step;
    return out;
}

double* shift(double* in, int size, int shftAmt) {
    // Shifts the numbers in "in" by"shftAmt" and the result is "out"
    double* out = new double[size];
    for (int i = 0; i < size; i++) {
        int ii = (i + shftAmt) % size;
        if (ii < 0) ii = ii + size;
        out[ii] = in[i];
    }
    return out;
}

double average(double* arr, int size) {
    double sum = 0;
    for (int i = 0; i < size; i++) 
        sum += arr[i];
    double avg = sum / size;
    return avg;
}

void setZero(double* arr, int size) {
    for (int i = 0; i < size; i++) 
        arr[i] = 0;
}
