function [S_vector_new, E_new] = ...
    function_Ising_rw_step_at_q(S_vector, N_atm, E, NN_table)
%
S_vector_new = S_vector;
%
% CHOOSE A RANDOM SPIN
flipped1_index = randsample(N_atm,1);
flipped1_pos_start = S_vector(flipped1_index,1);
flipped2_list2 = find( S_vector == - S_vector(flipped1_index,1));
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
% FLIP
%
S_vector_new(flipped1_index,1) = -flipped1_pos_start;
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
E_new = E - E_z_old_temp + E_z_new_temp;
%
% CHOOSE ANOTHER RANDOM SPIN THAT CAN BE FLIPPED TO KEEP STARTING M VALUE
%
% flipped2_list2 = find( SPM_new(1:N_atm) - SPM1_dif >= 1 & SPM_new(1:N_atm) - SPM1_dif <= Npos );
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
S_vector_new(flipped2_index,1) = -S_vector_new(flipped2_index,1);
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
end