function [S_vector, SPM, E, hist_nRPS, index_Npos_sum_vector] = ...
    function_WL_rw_criteria_correct(S_vector, SPM, E, S_vector_new, SPM_new, E_new, JDOS_nRPS_correct, E_list, hist_nRPS, index_Npos_sum_vector, index_Npos_sum_vector_new)
%
if rand < min( [ JDOS_nRPS_correct(E_list == E, index_Npos_sum_vector) ./ JDOS_nRPS_correct(E_list == E_new, index_Npos_sum_vector_new), 1]) % accept
    %
    E = E_new;
    S_vector = S_vector_new;
    SPM = SPM_new;
    index_Npos_sum_vector = index_Npos_sum_vector_new;
    %
    hist_nRPS(E_list == E_new, index_Npos_sum_vector_new) = hist_nRPS(E_list == E_new, index_Npos_sum_vector_new) + 1;
    %
else % reject
    %
    hist_nRPS(E_list == E, index_Npos_sum_vector) = hist_nRPS(E_list == E, index_Npos_sum_vector) + 1;
    %
end