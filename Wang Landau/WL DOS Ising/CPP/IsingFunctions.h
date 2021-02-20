#ifndef _ISING_FUNCTIONS_H
#define _ISING_FUNCTIONS_H

#include "IsingFunctions.cpp"

void randSpins(int*, int);
int CEnergy(int*, double*, double*, int, int, int);
int sumNeigh(int, int, int*, double*, double*, int, int);
void generateEnergies(int*, int);

#endif