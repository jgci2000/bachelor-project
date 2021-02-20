function [S_vector, E, hist_nRPS] = ...
    function_Ising_WL_rw_criteria_correct(S_vector, E, S_vector_new, E_new, JDOS_nRPS, E_list, hist_nRPS)
%
if rand < min( [ JDOS_nRPS(E_list == E, 1) ./ JDOS_nRPS(E_list == E_new, 1), 1]) % accept
    %
    E = E_new;
    S_vector = S_vector_new;
    %
    hist_nRPS(E_list == E_new, 1) = hist_nRPS(E_list == E_new, 1) + 1;
    %
else % reject
    %
    hist_nRPS(E_list == E, 1) = hist_nRPS(E_list == E, 1) + 1;
    %
end