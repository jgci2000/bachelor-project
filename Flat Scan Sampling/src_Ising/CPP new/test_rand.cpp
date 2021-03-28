#include <iostream>
#include "xoshiro256++.h"
#include "splimix64.h"
#include <chrono>
#include <stdint.h>

using namespace std;

int main()
{
	srand((unsigned) time(NULL));
	uint64_t seed = rand();
	splitmix64_seed(seed);

	for (int i = 0; i < 100; i++)
		cout << splitmix64() << " ";
	cout << endl;
		
	for (int i = 0; i < 4; i++)
		s[i] = splitmix64();

	for (int i = 0; i < 100; i++)
		cout << rand_xoshiro256pp() << " ";
	cout << endl;

	return 0;
}

