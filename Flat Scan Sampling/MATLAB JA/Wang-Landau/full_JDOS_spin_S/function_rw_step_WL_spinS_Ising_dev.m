function [E, M_z, S_vector, hist_WL, log_JDOS_WL] = ...
    function_rw_step_WL_spinS_Ising_dev(E, M_z, S_vector, hist_WL, N_atm, Npos, Z_spin_values, log_JDOS_WL, log_f, E_list, M_list, NN_table)
%
% RANDOM LINEAR INDEX OF FLIPPED SPIN
flip_pos = randi(N_atm,1);
%
% NEW RANDOM SPIN Z VALUE FOR FLIPPED SPIN
new_spin = Z_spin_values(randi(Npos,1));
%
% ENERGY CHANGE BY LOCAL UPDATE
%
% E_temp1 = -S_vector(flip_pos) .* (S_vector(NN_table(flip_pos,1)) + S_vector(NN_table(flip_pos,2)) + S_vector(NN_table(flip_pos,3)) + S_vector(NN_table(flip_pos,4)));
% E_temp2 = -new_spin .* (S_vector(NN_table(flip_pos,1)) + S_vector(NN_table(flip_pos,2)) + S_vector(NN_table(flip_pos,3)) + S_vector(NN_table(flip_pos,4)));
%

E_temp1 = 0;
E_temp2 = 0;
%
for a = 1:length(NN_table(1,:))
    %
    E_temp1 = E_temp1 - S_vector(flip_pos) .* (S_vector(NN_table(flip_pos,a)));
    E_temp2 = E_temp2 - new_spin .* (S_vector(NN_table(flip_pos,a)));
    %
end
%
E_new = E - E_temp1 + E_temp2;
%
% NEW MAGNETIZATION
M_z_new = M_z - S_vector(flip_pos) + new_spin;
%
% WL CRITERIA TO ACCEPT SPIN FLIP ATTEMPT
if rand < min( [ exp(log_JDOS_WL(E_list == E, M_list == M_z) - log_JDOS_WL(E_list == E_new, M_list == M_z_new)), 1]) % accept
    %
    E = E_new;
    M_z = M_z_new;
    S_vector(flip_pos) = new_spin;
    hist_WL(E_list == E, M_list == M_z) = hist_WL(E_list == E, M_list == M_z) + 1;
    log_JDOS_WL(E_list == E, M_list == M_z) = log_JDOS_WL(E_list == E, M_list == M_z) + log_f;
    %
else % reject
    %
    hist_WL(E_list == E, M_list == M_z) = hist_WL(E_list == E, M_list == M_z) + 1;
    log_JDOS_WL(E_list == E, M_list == M_z) = log_JDOS_WL(E_list == E, M_list == M_z) + log_f;
    %
end


%
end