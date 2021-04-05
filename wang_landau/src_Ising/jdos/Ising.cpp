//
// Source file for Ising spin 1/2 class.
// João Inácio, Mar. 25th, 2021
//
// Here are the functions definitions for the Ising header file.
// 


#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "Ising.h"
#include "WL_Functions.h"


void Ising::read_NN_talbe(std::string file_name)
{
    std::ifstream neighbour_tables_file(file_name);
    std::string line;

    if (neighbour_tables_file.is_open())
    {
        int i = 0;
        while (std::getline(neighbour_tables_file, line))
        {
            std::vector<std::string> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                this->NN_table[i++] = std::stold(a.at(idx));
        }
        neighbour_tables_file.close();
    }
    else
        std::cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << std::endl;
}


void Ising::read_norm_factor(std::string file_name)
{
    std::ifstream norm_factor_file(file_name);
    std::string line;

    if (norm_factor_file.is_open()) 
    {
        for (int i = 0; std::getline(norm_factor_file, line); i++)
            this->norm_factor[i] = std::stold(line);
        norm_factor_file.close();
    }
    else 
        std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
}

std::map<int, int> Ising::create_map(int init, int final, int step)
{
    std::map<int, int> out;
    int i = 0;
    while (init <= final)
    {
        out.insert(std::pair<int, int>(init, i));
        init += step;
        i++;
    }
    return out;
}

void Ising::set_E_config(int E_config)
{
    this->E_config = E_config;
    this->idx_E_config = this->energies[E_config];
}

void Ising::set_M_config(int M_config)
{
    this->M_config = M_config;
    this->idx_M_config = this->magnetizations[M_config];
}


std::vector<int> Ising::create_vector(int init, int final, int step)
{
    std::vector<int> out;
    while (init <= final)
    {
        out.push_back(init);
        init += step;
    }
    return out;
}


system_info Ising::get_system(int L, int lattice_num)
{
    system_info system;

    switch (lattice_num)
    {
        case 1:
            system.dim = 2;
            system.lattice = "SS";
            system.N_atm = L * L;
            system.NN = 4;
            break;

        case 2:
            system.dim = 3;
            system.lattice = "SC";
            system.N_atm = L * L * L;
            system.NN = 6;
            break;

        case 3:
            system.dim = 3;
            system.lattice = "BCC";
            system.N_atm = 2 * L * L * L; 
            system.NN = 8;
            break;
        
        case 4:
            system.dim = 3;
            system.lattice = "FCC";
            system.N_atm = 4 * L * L * L;
            system.NN = 12;
            break;

        case 5: 
            system.dim = 3;
            system.lattice = "HCP";
            system.N_atm = 2 * L * L * L;
            system.NN = 12;
            break;
        case 6:
            system.dim = 3;
            system.lattice = "Hex";
            system.N_atm = L * L * L;
            system.NN = 8;
            break;

        default:
            std::cout << "Invalid lattice number." << std::endl;
            break;
    }

    return system;
}




