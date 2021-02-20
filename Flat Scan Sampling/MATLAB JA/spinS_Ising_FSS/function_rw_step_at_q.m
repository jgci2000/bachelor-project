function [S_vector_new, SPM_new, E_new, index_Npos_sum_vector_new] = ...
    function_rw_step_at_q(S_vector, SPM, N_atm, Npos, E, NN_table, Z_spin_values, Npos_sum_configs)
%
% clear
% %
% S_vector = [2;2;2;2;2;2;2;2;2;2;2;2;2;0;2;2];
% SPM = [1;1;1;1;1;1;1;1;1;1;1;1;1;2;1;1];
% E = -64;
% Npos = 3;
% N_atm = 16;
% NN_table = [2,4,13,5;3,1,14,6;4,2,15,7;1,3,16,8;6,8,1,9;7,5,2,10;8,6,3,11;5,7,4,12;10,12,5,13;11,9,6,14;12,10,7,15;9,11,8,16;14,16,9,1;15,13,10,2;16,14,11,3;13,15,12,4];
% Z_spin_values = [2;0;-2];
%
S_vector_new = S_vector;
SPM_new = SPM;
%
% CHOOSE A RANDOM SPIN
flipped1_index = randsample(N_atm,1);
flipped1_pos_start = SPM(flipped1_index,1);
%
% CALC ENERGY BEFORE SPIN CHANGE
E_z_old_temp = 0;
%
for a = 1:length(NN_table(1,:))
    %
    E_z_old_temp = E_z_old_temp - S_vector_new(flipped1_index,1) .* S_vector_new(NN_table(flipped1_index,a),1);
    %
end
%
% NEW RANDOM, DIFFERENT SPIN VALUE
SPM1_start = SPM(flipped1_index,1);
SPM1_end_list = [1 : (flipped1_pos_start - 1), (flipped1_pos_start + 1) : Npos];
%
if length(SPM1_end_list) > 1
    %
    SPM1_end = randsample(SPM1_end_list,1);
    %
else
    %
    SPM1_end = SPM1_end_list;
    %
end
%
S_vector_new(flipped1_index,1) = Z_spin_values(SPM1_end);
SPM_new(flipped1_index,1) = SPM1_end;
%
% CALC NEW ENERGY
E_z_new_temp = 0;
%
for a = 1:length(NN_table(1,:))
    %
    E_z_new_temp = E_z_new_temp - S_vector_new(flipped1_index,1) .* S_vector_new(NN_table(flipped1_index,a),1);
    %
end
%
SPM1_dif = SPM1_end - SPM1_start;
E_new = E - E_z_old_temp + E_z_new_temp;
%
% CHOOSE ANOTHER RANDOM SPIN THAT CAN BE FLIPPED TO KEEP STARTING M VALUE
%
flipped2_list2 = find( SPM_new(1:N_atm) - SPM1_dif >= 1 & SPM_new(1:N_atm) - SPM1_dif <= Npos );
%
if length(flipped2_list2) > 1
    %
    flipped2_index = randsample(flipped2_list2,1);
    %
else
    %
    flipped2_index = flipped2_list2;
    %
end
%
% CALC ENERGY BEFORE SPIN CHANGE
E_z_old_temp = 0;
%
for a = 1:length(NN_table(1,:))
    %
    E_z_old_temp = E_z_old_temp - S_vector_new(flipped2_index,1) .* S_vector_new(NN_table(flipped2_index,a),1);
    %
end
%
% FLIP SPIN
SPM_new(flipped2_index) = SPM_new(flipped2_index) - SPM1_dif;
S_vector_new(flipped2_index,1) = Z_spin_values(SPM_new(flipped2_index));
%
% CALC NEW ENERGY
E_z_new_temp = 0;
%
for a = 1:length(NN_table(1,:))
    %
    E_z_new_temp = E_z_new_temp - S_vector_new(flipped2_index,1) .* S_vector_new(NN_table(flipped2_index,a),1);
    %
end
%
E_new = E_new - E_z_old_temp + E_z_new_temp;
%
index_Npos_sum_vector_new = find(sum(sum(SPM_new == 1:1:Npos) == Npos_sum_configs,2) == Npos);
%
end