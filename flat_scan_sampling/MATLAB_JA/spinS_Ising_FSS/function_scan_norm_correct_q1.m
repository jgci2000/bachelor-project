function [JDOS_nRPS_correct] = ...
    function_scan_norm_correct_q1(Npos, N_atm, E_list, NN, NN_table, Z_spin_values, Npos_sum_configs, JDOS_nRPS_correct)
%
SPM = ones(N_atm, 1);
%
%
prev_E_old_index = 1;
prev_Npos_sum_config_old_index = 1;
%
for x = 1:(Npos-1)
    %
    % FIND SPINS THAT CAN GO +x
    flip_list = find(SPM <= Npos -x);
    %
    % FLIP THEM SEQUENTIALLY, FIND NEW E
    for flip_index = 1:length(flip_list)
        %
        SPM_temp = SPM;
        %
        E_temp1 = 0;
        %
        for a = 1:length(NN_table(1,:))
            %
            E_temp1 = E_temp1 - Z_spin_values(SPM(flip_list(flip_index))) .* Z_spin_values(SPM(flip_list(flip_index)));%S_vector(NN_table(flip_list(flip_index),a));
            %
        end
        %
        E_temp2 = 0;
        %
        for a = 1:length(NN_table(1,:))
            %
            E_temp2 = E_temp2 - Z_spin_values(SPM(flip_list(flip_index)) + x) .* Z_spin_values(SPM(flip_list(flip_index)));%S_vector(NN_table(flip_list(flip_index),a));
            %
        end
        %
        % CALC E
        E_temp3 = max(Z_spin_values(:,1)).^2*(- N_atm * NN ./2) - E_temp1 + E_temp2;
        SPM_temp(flip_list(flip_index)) = SPM(flip_list(flip_index)) + x;
        %
        % UPDATE JDOS_correct
        index_Npos_sum_vector = find(sum(sum(SPM_temp == 1:1:Npos) == Npos_sum_configs{1+x,1},2) == Npos);
        JDOS_nRPS_correct{1+x,1}(E_list == E_temp3, index_Npos_sum_vector) = JDOS_nRPS_correct{1+x,1}(E_list == E_temp3, index_Npos_sum_vector) + JDOS_nRPS_correct{1,1}(prev_E_old_index,prev_Npos_sum_config_old_index);
        %
    end
    %
end
%
end
%
