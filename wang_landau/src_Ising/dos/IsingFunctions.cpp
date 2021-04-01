#include<iostream>
#include<time.h>


int sumNeigh(int i, int j, int*spins, double* plus, double* minus, int L, int nSpins) {
    // Computes the sum of the neighbours given (i, j)
    int pi = (int)plus[i] * L + j;
    if (pi == nSpins) pi = 0;
    int mi = (int)minus[i] * L + j;
    if (mi == nSpins) mi = 0;
    int pj = i * L + (int)plus[j];
    if (pj == nSpins) pj = 0;
    int mj = i * L + (int)minus[j];
    if (mj == nSpins) mj = 0;

    return spins[pi] + spins[mi] + spins[pj] + spins[mj];
}

int CEnergy(int* spins, double* plus, double* minus, int J, int L, int nSpins) {
    int E = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double S = (double) spins[i * L + j];
            double sumN = (double) sumNeigh(i, j, spins, plus, minus, L, nSpins);
            E += - J * sumN * S;
        }
    }
    return E / 2;
}

void generateEnergies(int* energies, int nSpins) {
    energies[0] = - 2 * nSpins;
    energies[nSpins - 2] = 2 * nSpins;
    for (int i = 1; i < nSpins - 2; i++) 
        energies[i] = 4 * (i + 1) - 2 * nSpins;
}

void randSpins(int* spins, int L) {
    // Gets a new seed each time we run
    srand((int)time(0));

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            // Generates a random integer (1 and 2)
            int rnd = (rand() % 2) + 1; 
            if (rnd == 1) spins[i * L + j] = -1;
            else spins[i * L + j] = +1;
        }
    }
}

