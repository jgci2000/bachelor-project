#ifndef _VECTOR_FUNCIIONS_H
#define _VECTOR_FUNCTIONS_H

#include "VectorFunctions.cpp"

void printSpins(int*);
void print1D(double*, int);
double* createVector(double, double, double);
double* shift(double*, int, int);
void setZero(double*, int);
double average(double*, int);

#endif