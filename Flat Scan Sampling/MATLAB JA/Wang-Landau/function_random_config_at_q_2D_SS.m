function [S_vector_WL, E_WL_old] = function_random_config_at_q_2D_SS(N_atm, NN, q, nnxpos, nnxneg, nnypos, nnyneg)
%
% INITIAL CONFIGURATION AT q
%
S_vector_WL = ones(N_atm, 1);
E_WL_old = -1 * N_atm * NN /2;
%
for index = 1:(q-1)
    %
    pos = find(S_vector_WL(:,1) == 1);
    flipped_pos = randsample(pos,1);
    %
    S_vector_WL(flipped_pos,1) = -1;
    %
    delta_E = - S_vector_WL(flipped_pos,1) .* ( ...
        S_vector_WL(nnxpos(flipped_pos),1) + ...
        S_vector_WL(nnxneg(flipped_pos),1) + ...
        S_vector_WL(nnypos(flipped_pos),1) + ...
        S_vector_WL(nnyneg(flipped_pos),1)) ; % energy of bonds to NN
    %
    E_WL_old = E_WL_old + 2*delta_E; % build the energy matrix
    %
end
