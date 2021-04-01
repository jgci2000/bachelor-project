clear
close all
%
% v44 - first version using WL rw from scratch, adapting v42
% v45 - cleanup, qol, more general stop condiftion for the WL rw
% v45 - *f or ^f WL DOS update choice and DOS keep choices
% v46 - qol
% v47 - no WL_DOS updates, skip enforced, equal to N_atm for Ising
% v48 - same as v47, cleanup, no use of external functions
%
% rng default
%
L = 4;
q_max = L^2/2 + 1;
%
max_rw_steps = 1E12;
skip = L^2; % 0 for no skip
REP = 1E4;
%
NN = 4;
N_atm = L^2;
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = -1/2 * N_atm * NN : +4 : 1/2 * N_atm * NN; % possible energy values
%
eval(['load norm_factor_Ising_Natm_', int2str(N_atm)])
%
% LATTICE AND NEIGHBOUR TABLES
%
lattice = repmat(1:L,L,1);
%
for index = 1:L
    %
    lattice(index,:) = lattice(index,:) + (index-1)*L;
    %
end
%
nnxpos = circshift(lattice,[0,-1]);
nnxneg = circshift(lattice,[0,1]);
nnypos = circshift(lattice,[1,0]);
nnyneg = circshift(lattice,[-1,0]);
%
nnxpos = reshape(nnxpos', [], 1);
nnxneg = reshape(nnxneg', [], 1);
nnypos = reshape(nnypos', [], 1);
nnyneg = reshape(nnyneg', [], 1);
%
JDOS_aprox = zeros(length(E_list), length(M_list));
JDOS_aprox_frac = zeros(length(E_list), length(M_list));
%
JDOS_aprox(1,1) = 1;
JDOS_aprox(1,length(M_list)) = 1;
JDOS_aprox_frac(1,1) = 1;
JDOS_aprox_frac(1,length(M_list)) = 1;
%
q = 1; % ONLY SCAN, USE ALL SPIN UP CONFIGS TO CALCULATE DOS(q=2)
%
neo_previous = zeros(length(E_list), length(E_list)); % old new
%
S_vector_old = ones(N_atm+1, 1);
E_old = -1/2 * N_atm * NN;
%
for k = 1:REP
    %
    S_vector_new = S_vector_old;
    pos = find(S_vector_old(:,1) == 1);
    %
    for spin_index = 1:N_atm
        %
        flipped_pos = pos(spin_index, 1);
        S_vector_new(flipped_pos,1) = -1; % flip the spin
        %
        delta_E = - S_vector_new(flipped_pos,1) .* ( ...
            S_vector_new(nnxpos(flipped_pos),1) + ...
            S_vector_new(nnxneg(flipped_pos),1) + ...
            S_vector_new(nnypos(flipped_pos),1) + ...
            S_vector_new(nnyneg(flipped_pos),1)) ; % energy of bonds to NN
        %
        E_new = E_old + 2*delta_E;
        %
        neo_previous(E_list == E_old, E_list == E_new) = neo_previous(E_list == E_old, E_list == E_new) + 1;
        %
        S_vector_new(flipped_pos,1) = +1; % unflip the spin
        %
    end
    %
end
%
% CALCULATE DOS at q = 2
%
for x = find(neo_previous)
    %
    [u,v] = ind2sub(length(E_list),x);
    JDOS_aprox_frac(v, q+1) = neo_previous(x)/sum(sum(neo_previous));
    %
end
%
hits = JDOS_aprox_frac(:, q) > 0;
%
% LOOP OVER q values, calculating DOS for q = q+1
%
for q = 2:(q_max-1)
    %
    neo_previous = zeros(length(E_list), length(E_list)); % old new
    %
    % WL RW USING neo_pad_dev_2_v6
    %
    WL_DOS = JDOS_aprox_frac(:, q); % DOS values start as the estimate from scan at q-1
    %
    WL_log_DOS = log(WL_DOS);
    %
    hist_WL = zeros(length(E_list),1);
    hist_E_selected_S_vector = zeros(length(E_list),1);
    %
    S_vector_WL = ones(N_atm, 1);
    %
    E_WL_old = -1/2 * N_atm * NN;
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
    hist_WL(E_list == E_WL_old) = hist_WL(E_list == E_WL_old) + 1;
    hist_E_selected_S_vector(E_list == E_WL_old) = hist_E_selected_S_vector(E_list == E_WL_old) + 1;
    %
    accept_counter = 1;
    reject_counter = 0;
    %
    hits = JDOS_aprox_frac(:, q) > 0;
    %
    k_saved = 1;
    k = 2;
    %
    while min(hist_E_selected_S_vector(hits)) < REP && min(hist_E_selected_S_vector(hist_E_selected_S_vector > 0)) < REP && k <= max_rw_steps
        %
        % choose an up and down spin to flip
        %
        pos = find(S_vector_WL(:,1) == 1);
        flipped_pos = randsample(pos,1);
        neg = find(S_vector_WL(:,1) == -1);
        %
        if length(neg) > 1
            flipped_neg = randsample(neg,1);
        else
            flipped_neg = neg;
        end
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
            if k >= k_saved + skip && hist_E_selected_S_vector(E_list == E_WL_new, 1) < REP
                %
                hist_E_selected_S_vector(E_list == E_WL_new, 1) = hist_E_selected_S_vector(E_list == E_WL_new, 1) + 1;
                k_saved = k;
                %
                % SCAN FOR WL ACCEPT
                %
                pos_scan = find(S_vector_WL(:,1) == 1);
                E_old = E_WL_old;
                %
                for spin_index = 1:length(pos_scan)
                    %
                    flipped_pos_scan = pos_scan(spin_index, 1);
                    S_vector_WL(flipped_pos_scan,1) = -1; % flip the spin
                    %
                    delta_E = - S_vector_WL(flipped_pos_scan,1) .* ( ...
                        S_vector_WL(nnxpos(flipped_pos_scan),1) + ...
                        S_vector_WL(nnxneg(flipped_pos_scan),1) + ...
                        S_vector_WL(nnypos(flipped_pos_scan),1) + ...
                        S_vector_WL(nnyneg(flipped_pos_scan),1)) ; % energy of bonds to NN
                    %
                    E_new = E_old + 2*delta_E;
                    %
                    neo_previous(E_list == E_old, E_list == E_new) = neo_previous(E_list == E_old, E_list == E_new) + 1;
                    %
                    S_vector_WL(flipped_pos_scan,1) = +1; % unflip the spin
                    %
                end
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
                pos_scan = find(S_vector_WL(:,1) == 1);
                E_old = E_WL_old;
                %
                for spin_index = 1:length(pos_scan)
                    %
                    flipped_pos_scan = pos_scan(spin_index, 1);
                    S_vector_WL(flipped_pos_scan,1) = -1; % flip the spin
                    %
                    delta_E = - S_vector_WL(flipped_pos_scan,1) .* ( ...
                        S_vector_WL(nnxpos(flipped_pos_scan),1) + ...
                        S_vector_WL(nnxneg(flipped_pos_scan),1) + ...
                        S_vector_WL(nnypos(flipped_pos_scan),1) + ...
                        S_vector_WL(nnyneg(flipped_pos_scan),1)) ; % energy of bonds to NN
                    %
                    E_new = E_old + 2*delta_E;
                    %
                    neo_previous(E_list == E_old, E_list == E_new) = neo_previous(E_list == E_old, E_list == E_new) + 1;
                    %
                    S_vector_WL(flipped_pos_scan,1) = +1; % unflip the spin
                    %
                end
                %
            end
            %
        end
        %
        k = k + 1;
        %
    end
    %
    % CALCULATE DOS at q = q+1
    %
    rehits = find(neo_previous);
    %
    for x = 1:length(rehits)
        %
        [u,v] = ind2sub(length(E_list),rehits(x));
        JDOS_aprox_frac(v, q+1) = JDOS_aprox_frac(v, q+1) + neo_previous(rehits(x))./sum(neo_previous(u,:)).*JDOS_aprox_frac(u,q);
        %
    end
    %
    JDOS_aprox(:,q) = JDOS_aprox_frac(:,q) .* norm_factor(q);
    %
end
%