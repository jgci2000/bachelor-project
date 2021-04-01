function [neo_previous, hist_E_selected_S_vector, hist_WL, accept_counter, reject_counter, k_saved, k, sample_ratio] = ...
    function_flat_test_inner3_v53(q, N_atm, NN, skip, E_list, REP, poolsize, JDOS_aprox_frac, max_rw_steps, nnxpos, nnxneg, nnypos, nnyneg)
%
neo_previous = zeros(length(E_list), length(E_list)); % old new
WL_log_DOS = log(JDOS_aprox_frac(:, q)); % DOS values start as the estimate from scan at q-1
%
hist_WL = zeros(length(E_list),1);
hist_E_selected_S_vector = zeros(length(E_list),1);
%
S_vector_WL = ones(N_atm, 1);
%
E_WL_old = -1/2 * N_atm * NN;
%
if q >= 2
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
    %
end
%
hist_WL(E_list == E_WL_old) = hist_WL(E_list == E_WL_old) + 1;
hist_E_selected_S_vector(E_list == E_WL_old) = hist_E_selected_S_vector(E_list == E_WL_old) + 1;
%
accept_counter = 1;
reject_counter = 0;
%
hits = JDOS_aprox_frac(:, q) > 0;
%
k_saved = 1;
%
% SCAN FOR FIRST CONFIGURATION
%
E_old = E_WL_old;
delta_E = S_vector_WL .* (S_vector_WL(nnxpos) + S_vector_WL(nnxneg) + S_vector_WL(nnypos) + S_vector_WL(nnyneg));
E_new = E_old + 2*delta_E(S_vector_WL(:,1) == 1);
neo_previous(E_list == E_old, :) = neo_previous(E_list == E_old, :) + sum(E_list == E_new',2)';
%
k = 2;
%
while min(hist_E_selected_S_vector(hits)) < floor(REP/poolsize) && min(hist_E_selected_S_vector(hist_E_selected_S_vector > 0)) < floor(REP/poolsize) && k <= max_rw_steps
    %
    % choose an up and down spin to flip
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
    E_WL_new = E_WL_old + 2*delta_E; % build the energy matrix
    %
    neg = find(S_vector_WL(:,1) == -1);
    %
    if length(neg) > 1
        flipped_neg = randsample(neg,1);
    else
        flipped_neg = neg;
    end
    %
    S_vector_WL(flipped_neg,1) = 1;
    %
    delta_E = - S_vector_WL(flipped_neg,1) .* ( ...
        S_vector_WL(nnxpos(flipped_neg),1) + ...
        S_vector_WL(nnxneg(flipped_neg),1) + ...
        S_vector_WL(nnypos(flipped_neg),1) + ...
        S_vector_WL(nnyneg(flipped_neg),1)) ; % energy of bonds to NN
    %
    E_WL_new = E_WL_new + 2*delta_E;
    %
    if rand < min( [ exp(WL_log_DOS(E_list == E_WL_old, 1) - WL_log_DOS(E_list == E_WL_new, 1)), 1]) % accept
        %
        E_WL_old = E_WL_new;
        hist_WL(E_list == E_WL_new, 1) = hist_WL( E_list == E_WL_new, 1) + 1; % update histogram
        %
        accept_counter = accept_counter + 1;
        %
        if k >= k_saved + skip && hist_E_selected_S_vector(E_list == E_WL_new, 1) < floor(REP/poolsize)
            %
            hist_E_selected_S_vector(E_list == E_WL_new, 1) = hist_E_selected_S_vector(E_list == E_WL_new, 1) + 1;
            k_saved = k;
            %
            % SCAN FOR WL ACCEPT
            %
            E_old = E_WL_old;
            delta_E = S_vector_WL .* (S_vector_WL(nnxpos) + S_vector_WL(nnxneg) + S_vector_WL(nnypos) + S_vector_WL(nnyneg));
            E_new = E_old + 2*delta_E(S_vector_WL(:,1) == 1);
            neo_previous(E_list == E_old, :) = neo_previous(E_list == E_old, :) + sum(E_list == E_new',2)';
            %
        end
        %
    else % reject
        %
        S_vector_WL(flipped_pos) = +1; % unflip spin
        S_vector_WL(flipped_neg) = -1; % unflip spin
        hist_WL(E_list == E_WL_old, 1) = hist_WL( E_list == E_WL_old, 1) + 1; % update histogram
        %
        reject_counter = reject_counter + 1;
        %
        if k >= k_saved + skip && hist_E_selected_S_vector(E_list == E_WL_old, 1) < REP
            %
            hist_E_selected_S_vector(E_list == E_WL_old, 1) = hist_E_selected_S_vector(E_list == E_WL_old, 1) + 1;
            k_saved = k;
            %
            % SCAN FOR WL REJECT
            %
            E_old = E_WL_old;
            delta_E = S_vector_WL .* (S_vector_WL(nnxpos) + S_vector_WL(nnxneg) + S_vector_WL(nnypos) + S_vector_WL(nnyneg));
            E_new = E_old + 2*delta_E(S_vector_WL(:,1) == 1);
            neo_previous(E_list == E_old, :) = neo_previous(E_list == E_old, :) + sum(E_list == E_new',2)';
            %
        end
        %
    end
    %
    k = k + 1;
    %
end
%
sampled = REP/poolsize * nnz(hits);
sample_ratio = sampled/k * skip;
%
end
%