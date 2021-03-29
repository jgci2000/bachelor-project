//
// Header file for Ising spin 1/2 class.
// João Inácio, Mar. 25th, 2021
//
// Here are the class structure, atributes and functions declarations.
//


#ifndef ISING_H
#define ISING_H


struct system_info 
{
    int N_atm;
    int NN;
    std::string lattice;
    int dim;
};


class Ising {
    private:
        int lattice_num;

        std::vector<int> create_vector(int, int, int);
        std::map<int, int> create_map(int, int, int);
        system_info get_system(int, int);

    public:
        int L;
        int N_atm;
        int NN;
        int J = 1;

        int NE;
        int NM;
        int max_E;
        int max_M;
        
        std::map<int, int> energies;
        std::map<int, int> magnetizations;

        int dim;
        std::string lattice;

        int *spins_vector;
        int *NN_table;
        long double *norm_factor;

        int E_config;
        int idx_E_config;

        int M_config;
        int idx_M_config;

        Ising(int L, int lattice_num)
        {
            system_info system = get_system(L, lattice_num);

            this->L = L;
            this->lattice = system.lattice;
            this->dim = system.dim;
            this->NN = system.NN;
            this->N_atm = system.N_atm;

            this->max_E = (1.0 / 2.0) * this->NN * this->N_atm;
            this->max_M = this->N_atm;

            this->NE = 1 + (this->max_E / 2);
            this->NM = this->N_atm + 1;

            this->spins_vector = new int[this->N_atm];
            this->NN_table = new int[this->N_atm * this->NN];
            this->norm_factor = new long double[this->NM];

            this->energies = create_map(- max_E, max_E, 4);
            this->magnetizations = create_map(- max_M, max_M, 2);
        }

        void read_NN_talbe(std::string);
        void read_norm_factor(std::string);

        void set_E_config(int);
        void set_M_config(int);

        ~Ising()
        {
            delete[] this->spins_vector, this->NN_table, this->norm_factor;
        }

};

#endif // ISING_H
