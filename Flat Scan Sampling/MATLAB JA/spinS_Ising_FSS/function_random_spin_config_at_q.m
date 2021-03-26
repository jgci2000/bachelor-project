function [S_vector, SPM, E, index_Npos_sum_vector] = ...
    function_random_spin_config_at_q(N_atm, NN, Z_spin_values, q, Npos, NN_table, Npos_sum_configs)
%
S_vector = repmat(Z_spin_values(1,1), [N_atm 1]);
SPM = ones(N_atm,1);
%
E = -1/2 * NN * N_atm * Z_spin_values(1,1).^2;
%
for q_function = 2:q
    %
    flipped_pos = find(SPM < Npos);
    flipped_pos_index = randsample(flipped_pos,1);
    %flipped_pos_index = randi(length(flipped_pos));
    %
    % CALC ENERGY BEFORE SPIN CHANGE
    E_z_old_temp = 0;
    %
    for a = 1:length(NN_table(1,:))
        %
        E_z_old_temp = E_z_old_temp - S_vector(flipped_pos_index,1) .* S_vector(NN_table(flipped_pos_index,a),1);
        %
    end
    %
    % FLIP SPIN
    SPM(flipped_pos_index,1) = SPM(flipped_pos_index,1) + 1;
    S_vector(flipped_pos_index,1) = Z_spin_values(SPM(flipped_pos_index,1),1);
    %
    % CALC NEW ENERGY
    E_z_new_temp = 0;
    %
    for a = 1:length(NN_table(1,:))
        %
        E_z_new_temp = E_z_new_temp - S_vector(flipped_pos_index,1) .* S_vector(NN_table(flipped_pos_index,a),1);
        %
    end
    %
    E = E - E_z_old_temp + E_z_new_temp;
    %
end
%
index_Npos_sum_vector = find(sum(sum(SPM == 1:1:Npos) == Npos_sum_configs,2) == Npos);
% index_Npos_sum_vector = find(sum(sum(SPM_temp == 1:1:Npos) == Npos_sum_configs{1+x,1},2) == Npos);
%
end