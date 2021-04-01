/*
    Wang Landau Sampling for the 2D Ising Model
    João Inácio, Aug. 22, 2020
*/
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <algorithm>

void printSpins(int*, int);
void print1D(double*, int);
void randSpins(int*);
int CEnergy(int*, double*, double*);
int sumNeigh(int, int, int*, double*, double*);
double* createVector(double, double, double);
double* shift(double*, int, int);
int* generateEnergies();
double average(double*, int);
void normalizeDOS(int*, double*, int);
void setZero(double*, int);

const int J = 1;                // Interaction strength
const int L = 2;                // Lattice size
const int N_SPINS = L * L;      // Number of particles

int main() {
    // Start measuring time
    clock_t start = clock();

    double* aux = createVector(0, L-1, 1);
    double* minus = shift(aux, L, -1);
    double* plus = shift(aux, L, 1);

    double f = exp(1);          // Modification factor
    double p = 0.95;            // Flatness -> min(Hist) > avg(Hist)*p

    // All of the possible configurations
    int* energies = generateEnergies();
    int NE = N_SPINS - 1;
    
    // Create spins on the heap so it can be modified in the funcions
    int* spins = new int[L * L];
    randSpins(spins);
    // printSpins(spins, L);
    int E = CEnergy(spins, plus, minus); // Energy of the first configuration
    auto itr = std::find(energies, energies + NE, E);
    int idxE = std::distance(energies, itr);

    // DOS array and Histogram
    double* lngE = new double[NE];
    setZero(lngE, NE);
    double* hist = new double[NE];
    setZero(hist, NE);
    
    // MCSweeps counter
    int MCSweeps = 0;

    // Random seed
    srand((int)time(0));

    while(f > (1 + pow(10, -8))) {
        // Each loop is one MCSweep
        for (int ni = 0; ni < N_SPINS; ni++) {
            int i = (rand() % L); // Generates a random int (1 - L)
            int j = (rand() % L);

            // dE = 2 * J * S * sum(neighbours) and dE = ENew - E, so
            int S = spins[i + j * L];
            int sumN = sumNeigh(i, j, spins, plus, minus);
            int ENew = E + 2 * J * S * sumN;

            auto itr = std::find(energies, energies + NE, ENew);
            int idxENew = std::distance(energies, itr);

            // Flip the spin
            double ratio = exp(lngE[idxE] - lngE[idxENew]);
            double P = std::min(ratio, 1.0);

            if (P > ((double) rand() / (RAND_MAX))) {
                spins[i + j * L] = - S;
                E = ENew;
                idxE = idxENew;
            }
            // printSpins(spins);
            // Update the histogram and g(E)
            hist[idxE] += 1;
            lngE[idxE] += log(f);
            // print1D(lngE, NE);
        }

        MCSweeps += 1;

        // Check at each 10000 sweeps
        if (MCSweeps % 10000 == 0) {
            double avgH = average(hist, NE);
            int minH = *std::min_element(hist, hist + NE);
            // std::cout << MCSweeps << std::endl;
            // for (int i = 0; i < NE; i++) 
            //     std::cout << hist[i] << " : " << energies[i] << std::endl;
            // std::cout << avgH << std::endl;
            if (minH > avgH * p) {
                printf("%d: the histogram is flat. Min: %d Avg: %.1f f: %.8f\n", MCSweeps, minH, avgH, f);
                f = sqrt(f);
                MCSweeps = 0;
            }
        }
    }

    // Normalize the DOS
    normalizeDOS(energies, lngE, NE);
    print1D(lngE, NE);

    // Erase arrays from the heap
    delete[] spins, aux, minus, plus, energies, hist, lngE;

    // Stop measuring time and calculate the elapsed time
    clock_t end = clock();
    double elapsed = double(end - start)/CLOCKS_PER_SEC;
    printf("Time measured: %.7f seconds.\n", elapsed);
    return 0;
}

int CEnergy(int* spins, double* plus, double* minus) {
    int E = 0;
    for (int j = 0; j < L; j++) {
        for (int i = 0; i < L; i++) {
            double S = (double) spins[i + j * L];
            double sumN = sumNeigh(i, j, spins, plus, minus);
            E += - J * sumN * S;
        }
    }
    return E / 2;
}

void normalizeDOS(int* energies, double* lngE, int NE) {
    auto itr = std::find(energies, energies + NE, - 2  * N_SPINS);
    int idx = std::distance(energies, itr);

    for (int i = 0; i < NE; i++) 
        lngE[i] = lngE[i] - lngE[idx] + log(2);
}

int sumNeigh(int i, int j, int*spins, double* plus, double* minus) {
    // Computes the sum of the neighbours given (i, j)
    int pi = (int)plus[i] + j * L;
    if (pi == N_SPINS) pi = 0;
    int mi = (int)minus[i] + j * L;
    if (mi == N_SPINS) mi = 0;
    int pj = i + (int)plus[j] * L;
    if (pj == N_SPINS) pj = 0;
    int mj = i + (int)minus[j] * L;
    if (mj == N_SPINS) mj = 0;

    return spins[pi] + spins[mi] + spins[pj] + spins[mj];
}

int* generateEnergies() {
    int* energies = new int[N_SPINS - 1];

    energies[0] = - 2 * N_SPINS;
    energies[N_SPINS - 2] = 2 * N_SPINS;
    for (int i = 1; i < N_SPINS - 2; i++) 
        energies[i] = 4 * (i + 1) - 2 * N_SPINS;

    return energies;
}

void randSpins(int* spins) {
    // Gets a new seed each time we run
    srand((int)time(0));

    for (int j = 0; j < L; j++) {
        for (int i = 0; i < L; i++) {
            // Generates a random integer (1 and 2)
            int rnd = (rand() % 2) + 1; 
            if (rnd == 1) spins[i + j * L] = -1;
            else spins[i + j * L] = +1;
        }
    }
}

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
