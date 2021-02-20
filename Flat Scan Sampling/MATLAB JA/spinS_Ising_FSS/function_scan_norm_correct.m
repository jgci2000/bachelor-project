function [JDOS_nRPS_correct] = ...
    function_scan_norm_correct(S_vector, SPM, E, Npos, E_list, NN_table, Z_spin_values, REP, Npos_sum_configs, JDOS_nRPS_correct, q)
%
prev_E_old_index = E_list == E;
prev_Npos_sum_config_old_index = sum(sum(SPM == 1:1:Npos) == Npos_sum_configs{q,:},2) == Npos;
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
            E_temp1 = E_temp1 - S_vector(flip_list(flip_index),1) .* S_vector(NN_table(flip_list(flip_index),a));
            %
        end
        %
        E_temp2 = 0;
        %
        for a = 1:length(NN_table(1,:))
            %
            E_temp2 = E_temp2 - Z_spin_values(SPM(flip_list(flip_index)) + x) .* S_vector(NN_table(flip_list(flip_index),a));
            %
        end
        %
        % CALC E
        E_temp3 = E - E_temp1 + E_temp2;
        SPM_temp(flip_list(flip_index)) = SPM(flip_list(flip_index)) + x;
        %
        % UPDATE JDOS_correct
        index_Npos_sum_vector = find(sum(sum(SPM_temp == 1:1:Npos) == Npos_sum_configs{q+x,1},2) == Npos);
        JDOS_nRPS_correct{q+x,1}(E_list == E_temp3, index_Npos_sum_vector) = JDOS_nRPS_correct{q+x,1}(E_list == E_temp3, index_Npos_sum_vector) + JDOS_nRPS_correct{q,1}(prev_E_old_index, prev_Npos_sum_config_old_index)/REP;
        %
    end
    %
end
%
end