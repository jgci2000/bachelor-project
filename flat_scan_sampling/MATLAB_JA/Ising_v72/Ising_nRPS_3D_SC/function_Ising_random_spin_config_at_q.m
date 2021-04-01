function [S_vector, E] = ...
    function_Ising_random_spin_config_at_q(N_atm, NN, q, NN_table)
%
S_vector = ones(N_atm,1);
%
E = - NN * N_atm /2; 
%
for q_function = 2:q
    %
    flipped_pos = find(S_vector == 1);
    flipped_pos_index = randsample(flipped_pos,1);
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
    S_vector(flipped_pos_index,1) = -1;
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
end