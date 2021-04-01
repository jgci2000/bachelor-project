clear
close all
%
% v44 - first version using WL rw from scratch, adapting v42
% v45 - cleanup, qol, more general stop condiftion for the WL rw
% v45 - *f or ^f WL DOS update choice and DOS keep choices
% v46 - qol
% v47 - no WL_DOS updates, skip enforced, equal to N_atm for Ising
% v48 - same as v47, cleanup, no use of external functions
% v49 - vectorized scan - no real speedup
% v50 - shorten code
%
rng default
%
L = 4;
q_max = L^2 / 2 + 1;
%
max_rw_steps = 1E12;
skip = L^2; % 0 for no skip
REP = 1E3;
%
NN = 4;
N_atm = L^2;
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = -1/2 * N_atm * NN : +4 : 1/2 * N_atm * NN; % possible energy values
%
wspace_filename = ['workspace_ft_v50_L', int2str(L), 'R1E', int2str(log10(REP))];
%
if q_max == L^2 + 1
    %
    wspace_filename = [wspace_filename, '_q_full'];
    disp('full sweeps from M to -M')
    %
end
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', wspace_filename]);
disp(['output workspace filename: ', wspace_filename])
%
disp(['L = ', int2str(L)])
disp(['REP = 1E', num2str(log10(REP))])
disp(['skip = ', int2str(skip)])
%
eval(['load norm_factor_Ising_Natm_', int2str(N_atm)])
%
% compare to exact JDOS of 4x4 system
%
if L == 4
    %
    load JDOS_L4_2D_SS_exact.mat
    load norm_factor_Ising_Natm_16.mat
    %
    JDOS_exact_frac = JDOS_exact;
    %
    for q_x = 1:(N_atm+1)
        %
        JDOS_exact_frac(:,q_x) = JDOS_exact(:,q_x) ./ sum(JDOS_exact(:,q_x));
        %
    end
    %
    JDOS_erro = zeros(length(E_list), q_max);
    JDOS_erro_total = 0;
    JDOS_erro_max = zeros(q_max, 1);
    %
    output = nan(q_max,7);
    %
else
    %
    output = nan(q_max, 5);
    %
end
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
% LOOP OVER q values, calculating DOS for q = q+1
%
tStart = tic;
%
for q = 1:(q_max-1)
    %
    tic
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
    k = 2;
    %
    while min(hist_E_selected_S_vector(hits)) < REP && min(hist_E_selected_S_vector(hist_E_selected_S_vector > 0)) < REP && k <= max_rw_steps
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
            if k >= k_saved + skip && hist_E_selected_S_vector(E_list == E_WL_new, 1) < REP
                %
                hist_E_selected_S_vector(E_list == E_WL_new, 1) = hist_E_selected_S_vector(E_list == E_WL_new, 1) + 1;
                k_saved = k;
                %
                % SCAN FOR WL ACCEPT
                %
                pos_scan = find(S_vector_WL(:,1) == 1);
                E_old = E_WL_old;
                delta_E = S_vector_WL .* (S_vector_WL(nnxpos) + S_vector_WL(nnxneg) + S_vector_WL(nnypos) + S_vector_WL(nnyneg));
                E_new = E_old + 2*delta_E(pos_scan);               
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
                pos_scan = find(S_vector_WL(:,1) == 1);
                E_old = E_WL_old;
                delta_E = S_vector_WL .* (S_vector_WL(nnxpos) + S_vector_WL(nnxneg) + S_vector_WL(nnypos) + S_vector_WL(nnyneg));
                E_new = E_old + 2*delta_E(pos_scan);
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
    JDOS_aprox(:,q+1) = JDOS_aprox_frac(:,q+1) .* norm_factor(q+1);
    %
    q_timer = toc;
    %
    if L ~= 4
        %
        disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', num2str(nnz(hits)), ' | time per E pt: ', num2str(q_timer/nnz(hits)), ' secs'])
        %
        output(q, 1:5) = [q, q_max, q_timer, nnz(hits), q_timer/nnz(hits)];
        %
    else
        %
        JDOS_erro(JDOS_exact_frac(:,q) ~= 0, q) = (JDOS_aprox_frac(JDOS_exact_frac(:,q) ~= 0,q) - JDOS_exact_frac(JDOS_exact_frac(:,q) ~= 0, q)) ./ JDOS_exact_frac(JDOS_exact_frac(:,q) ~= 0, q);
        JDOS_erro_total = JDOS_erro_total + sum(abs(JDOS_erro(:,q)));
        JDOS_erro_max(q,1) = max(abs(JDOS_erro(:,q)));
        %
        disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', num2str(nnz(hits)), ' | time per E pt: ', num2str(q_timer/nnz(hits)), ' secs | JDOS_err_q: ', num2str(sum(abs(JDOS_erro(:,q)))), ' | cum JDOS_err: ', num2str(JDOS_erro_total)]);
        %
        output(q, 1:7) = [q, q_max, q_timer, nnz(hits), q_timer/nnz(hits), sum(abs(JDOS_erro(:,q))), JDOS_erro_total];
        %
    end
    %
end
%
tEnd = toc(tStart);
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(tEnd), ' secs']);
%
if q_max == (L^2/2 + 1) || q_max == N_atm + 1
    %
    eval(['save ', wspace_filename, ' -v7.3'])
    %
end
