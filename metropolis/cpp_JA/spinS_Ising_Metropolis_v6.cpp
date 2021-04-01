// v1 - works for multiple T and H values, checked with Matlab output
// v2 - output data as individual ascii files
// v3 - speedup (no vectors)
// v4 - define function
// v5 - QoL output run_info.txt
// v6 - corrections for non-integer M values

#include <stdint.h>
#include <armadillo>
#include <ctime>
#include "mt19937ar.h"

using namespace std;
using namespace arma;

rowvec myFunction(double T_temp, double H_temp, int N_atm, int NN, int Npos, colvec Z_spin_values, imat NN_table, int SHF, int REP_interval, int REP_samples, int skip) {
  
    // cout << T_temp << " " << H_temp << endl;

    long double M_for_averaging = 0;
    long double M2_for_averaging = 0;
    long double M4_for_averaging = 0;
    long double E_for_averaging = 0;
    long double E2_for_averaging = 0;
    long double abs_M_for_averaging = 0;

    // ALL SPINS UP

    colvec S_vector(N_atm, fill::ones);
    double Mag = 1.0 * N_atm;
    double E = 1.0 * (-N_atm * NN / 2 -Mag*H_temp);

    // SHUFFLE

    colvec S_vector_new(N_atm, fill::ones);

    int flipped1_index;
    int k_shf;
    for (k_shf = 0; k_shf < SHF; k_shf++) {

        // TRIAL SPIN CHANGE

        S_vector_new = S_vector;
        flipped1_index = floor(genrand_res53()*N_atm);
        double E_z_old_temp = 0;
        int a;

        for (a = 0; a < NN; a++)
            E_z_old_temp = E_z_old_temp - S_vector_new(flipped1_index) * S_vector_new(NN_table(flipped1_index,a)-1);
        E_z_old_temp = E_z_old_temp - Mag*H_temp;

        S_vector_new(flipped1_index) = Z_spin_values(floor(genrand_res53()*Npos));
        double M_new = Mag - S_vector(flipped1_index) + S_vector_new(flipped1_index);
        double E_z_new_temp = 0;

        for (a = 0; a < NN; a++)
            E_z_new_temp = E_z_new_temp - S_vector_new(flipped1_index) * S_vector_new(NN_table(flipped1_index,a)-1);
        E_z_new_temp = E_z_new_temp - M_new*H_temp;

        double E_new = E - E_z_old_temp + E_z_new_temp;

        // METROPOLIS RW CRITERIA
        
        if (genrand_res53() < exp(-(E_new - E) / T_temp)) {
            E = E_new;
            S_vector = S_vector_new;
            Mag = M_new;
        }
    }

    // TRIAL FLIP, ACCUMULATE FOR AVERAGING

    int k_rep;

    for (k_rep = 0; k_rep <= REP_interval; k_rep++) {

        // TRIAL SPIN CHANGE

        S_vector_new = S_vector;
        flipped1_index = floor(genrand_res53()*N_atm);
        double E_z_old_temp = 0;
        int a;

        for (a = 0; a < NN; a++)
            E_z_old_temp = E_z_old_temp - S_vector_new(flipped1_index) * S_vector_new(NN_table(flipped1_index,a)-1);
        E_z_old_temp = 1.0 * (E_z_old_temp - Mag*H_temp);

        S_vector_new(flipped1_index) = Z_spin_values(floor(genrand_res53()*Npos));
        double M_new = Mag - S_vector(flipped1_index) + S_vector_new(flipped1_index);
        double E_z_new_temp = 0;

        for (a = 0; a < NN; a++)
            E_z_new_temp = E_z_new_temp - S_vector_new(flipped1_index) * S_vector_new(NN_table(flipped1_index,a)-1);

        E_z_new_temp = E_z_new_temp - M_new*H_temp;
        double E_new = E - E_z_old_temp + E_z_new_temp;

        // METROPOLIS RW CRITERIA

        if (genrand_res53() < exp(-(E_new - E) / T_temp)) {
            E = E_new;
            S_vector = S_vector_new;
            Mag = M_new;
        }

        if (k_rep > 0 && k_rep % skip == 0) {

            M_for_averaging = M_for_averaging + Mag;
            M2_for_averaging = M2_for_averaging + pow(Mag,2);
            M4_for_averaging = M4_for_averaging + pow(Mag,4);
            E_for_averaging = E_for_averaging + E;
            E2_for_averaging = E2_for_averaging + pow(E,2);
            abs_M_for_averaging = abs_M_for_averaging + abs(Mag);
        }

    }

    double avg_M = 1.0 * M_for_averaging / REP_samples;
    double avg_M2 = 1.0 * M2_for_averaging / REP_samples;
    double avg_M4 = 1.0 * M4_for_averaging / REP_samples;
    double avg_E = E_for_averaging / REP_samples;
    double avg_E2 = E2_for_averaging / REP_samples;
    double avg_abs_M = 1.0 * abs_M_for_averaging / REP_samples;

    // cout << avg_M << " " << avg_M2 << " " <<  avg_M4 << " " <<  avg_E << " " << avg_E2 << " " << avg_abs_M << endl;

    return {avg_M, avg_M2, avg_M4, avg_E, avg_E2, avg_abs_M};

}

int main(void)
{

    uint32_t seed = 1;
    init_genrand(seed);

    time_t tstart, tend; 
    tstart = time(0);

    colvec T(8, fill::none);
    int i;
    for (i = 0; i < 8; i++)
        T[i] = 1.0 * (8 - i*2);

    colvec H(11, fill::none);
    int j;
    for (j = 0; j < 11; j++)
        H[j] = 1.0 * (1 - j*0.1);

    int L = 16;
    int Npos = 2; // (2*S)+1

    int SHF = 1E5;
    int REP_interval = 1E5;
    int REP_samples = 1E5;

    int N_atm = pow(L,3);
    int NN = 8;

    int T_length = T.n_elem;
    int H_length = H.n_elem;

    double spin = 1.0 * (Npos-1)/2;

    if (Npos > 2)
    	T = T / ((0.5+1)/0.5) * ((spin+1)/spin); // auto re-scale of T for S > 1/2

    int skip = REP_interval/REP_samples;

    imat NN_table(N_atm, NN, fill::none);
    NN_table.load("neighbour_table_3D_Hex_1NN_L16.txt", raw_ascii);

    colvec Z_spin_values(Npos, fill::ones);
    for (i = 0; i <= Npos; ++i) 
        Z_spin_values[i] = 1.0 * (-(Npos-1) + i*2);
    Z_spin_values = 1.0 * Z_spin_values / (Npos-1);

    cube output(T_length, 9, H_length, fill::none);

    for (j = 0; j < H_length; j++) {

        for (i = 0; i < T_length; i++) {

            output(span(i,i), span(0,5), span(j,j)) = myFunction(T[i], H[j], N_atm, NN, Npos, Z_spin_values, NN_table, SHF, REP_interval, REP_samples, skip); // call the function

            cout << "T: " << T[i] << endl;

        }

        output.slice(j).col(6) = (output.slice(j).col(4) - pow(output.slice(j).col(3),2)) / (pow(T,2)); // Cp
        output.slice(j).col(7) = (output.slice(j).col(1) - pow(output.slice(j).col(0),2)) / T; // susc
        output.slice(j).col(8) = 1 - output.slice(j).col(2) / (3 * pow(output.slice(j).col(1),2)); // BinderC

        // cout << output.slice(j) << endl;

        cout << "H: " << H[j] << endl;

    }

    // cout << output.col(0) << endl;

    mat avg_M_mat_temp(H_length, T_length, fill::none);
    avg_M_mat_temp = reshape(output.col(0), T_length, H_length, 1);
    mat avg_M_mat(H_length, T_length, fill::none);
    avg_M_mat = avg_M_mat_temp.t();
    // cout << avg_M_mat << endl;
    avg_M_mat.save("avg_M.dat", raw_ascii);

    mat avg_M2_mat_temp(H_length, T_length, fill::none);
    avg_M2_mat_temp = reshape(output.col(1), T_length, H_length, 1);
    mat avg_M2_mat(H_length, T_length, fill::none);
    avg_M2_mat = avg_M2_mat_temp.t();
    // cout << avg_M2_mat << endl;
    avg_M2_mat.save("avg_M2.dat", raw_ascii);

    mat avg_M4_mat_temp(H_length, T_length, fill::none);
    avg_M4_mat_temp = reshape(output.col(2), T_length, H_length, 1);
    mat avg_M4_mat(H_length, T_length, fill::none);
    avg_M4_mat = avg_M4_mat_temp.t();
    // cout << avg_M4_mat << endl;
    avg_M4_mat.save("avg_M4.dat", raw_ascii);

    mat avg_E_mat_temp(H_length, T_length, fill::none);
    avg_E_mat_temp = reshape(output.col(3), T_length, H_length, 1);
    mat avg_E_mat(H_length, T_length, fill::none);
    avg_E_mat = avg_E_mat_temp.t();
    // cout << avg_E_mat << endl;
    avg_E_mat.save("avg_E.dat", raw_ascii);

    mat avg_E2_mat_temp(H_length, T_length, fill::none);
    avg_E2_mat_temp = reshape(output.col(4), T_length, H_length, 1);
    mat avg_E2_mat(H_length, T_length, fill::none);
    avg_E2_mat = avg_E2_mat_temp.t();
    // cout << avg_E2_mat << endl;
    avg_E2_mat.save("avg_E2.dat", raw_ascii);

    mat avg_abs_M_mat_temp(H_length, T_length, fill::none);
    avg_abs_M_mat_temp = reshape(output.col(5), T_length, H_length, 1);
    mat avg_abs_M_mat(H_length, T_length, fill::none);
    avg_abs_M_mat = avg_abs_M_mat_temp.t();
    // cout << avg_abs_M_mat << endl;
    avg_abs_M_mat.save("avg_abs_M.dat", raw_ascii);

    mat avg_Cp_mat_temp(H_length, T_length, fill::none);
    avg_Cp_mat_temp = reshape(output.col(6), T_length, H_length, 1);
    mat avg_Cp_mat(H_length, T_length, fill::none);
    avg_Cp_mat = avg_Cp_mat_temp.t();
    // cout << avg_Cp_mat << endl;
    avg_Cp_mat.save("avg_Cp.dat", raw_ascii);

    mat avg_susc_mat_temp(H_length, T_length, fill::none);
    avg_susc_mat_temp = reshape(output.col(7), T_length, H_length, 1);
    mat avg_susc_mat(H_length, T_length, fill::none);
    avg_susc_mat = avg_susc_mat_temp.t();
    // cout << avg_susc_mat << endl;
    avg_susc_mat.save("avg_susc.dat", raw_ascii);

    mat avg_BinderC_mat_temp(H_length, T_length, fill::none);
    avg_BinderC_mat_temp = reshape(output.col(8), T_length, H_length, 1);
    mat avg_BinderC_mat(H_length, T_length, fill::none);
    avg_BinderC_mat = avg_BinderC_mat_temp.t();
    // cout << avg_BinderC_mat << endl;
    avg_BinderC_mat.save("avg_BinderC.dat", raw_ascii);

    T.save("T.dat", raw_ascii);
    H.save("H.dat", raw_ascii);

    tend = time(0); 
    cout << "calc time:" << difftime(tend, tstart) << " second(s)." << endl;

    ofstream cout("run_info.txt");

    cout << "L: " << L << endl;
    cout << "Npos: " << Npos << endl;
    cout << "N_atm: " << N_atm << endl;
    cout << "NN: " << NN << endl;
    cout << "SHF: " << SHF << " = 10^" << log10(SHF) << endl;
    cout << "REP_interval: " << REP_interval << " = 10^" << log10(REP_interval) << endl;
    cout << "REP_samples: " << REP_samples << " = 10^" << log10(REP_samples) << endl;
    cout << "T_length: " << T_length << endl;
    cout << "H_length: " << H_length << endl;
    cout << "calc time:" << difftime(tend, tstart) << " second(s)." << endl;

    return 0;
}