#include <iostream>
#include <array>
#include <vector>

#include "Ising.h"

using namespace std;

void do_something(Ising &ising, array<vector<int>, 2> &flip_list)
{
    flip_list[0].pop_back();
    flip_list[1].push_back(15);

    ising.spins_vector[ising.N_atm - 1] = - 1;
    
    cout << "inside function: " << endl;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < flip_list[i].size(); j++)
            cout << flip_list[i].at(j) << " ";
        cout << endl;
    }

    for (int i = 0; i < ising.N_atm; i++)
        cout << ising.spins_vector[i] << " ";
    cout << endl;
}

int main()
{
    Ising ising(4, 1);

    for (int i = 0; i < ising.N_atm; i++)
        ising.spins_vector[i] = 1;
    ising.set_E_config(- ising.max_E);
    
    array<vector<int>, 2> flip_list;
    for (int i = 0; i < ising.N_atm; i++)
        flip_list[0].push_back(i);
    
    do_something(ising, flip_list);

    cout << "outside function: " << endl;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < flip_list[i].size(); j++)
            cout << flip_list[i].at(j) << " ";
        cout << endl;
    }

    for (int i = 0; i < ising.N_atm; i++)
        cout << ising.spins_vector[i] << " ";
    cout << endl;

    return 0;
}

