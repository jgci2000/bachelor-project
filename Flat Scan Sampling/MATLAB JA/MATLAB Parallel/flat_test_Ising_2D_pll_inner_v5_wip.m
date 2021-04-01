clear
close all
%
% v1 - from v38
% v2 - cleanup
% v3 - full magnetization (positive and negative) sampling
% v4 - update from v50
% v5 - update from v53; more outputs, add core# to output filename
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
p = gcp('nocreate'); % If no pool, do not create new one.
%
if isempty(p)
    poolsize = 1;
else
    poolsize = p.NumWorkers;
end
%
NN = 6;
N_atm = L^2;
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = -1/2 * N_atm * NN : +4 : 1/2 * N_atm * NN; % possible energy values
%
wspace_filename = ['workspace_ft_Ising_2D_pll_inner_v5_L', int2str(L), 'R1E', int2str(log10(REP))];
%
if skip ~= N_atm
    %
    wspace_filename = [wspace_filename, '_skip', int2str(skip)];
    disp('skip parameter not equal to number of spins')
    %
end
%
if q_max == N_atm
    %
    q_max = N_atm - 1;
    %
end
%
if q_max == N_atm - 1
    %
    wspace_filename = [wspace_filename, '_q_full'];
    disp('full sweeps from M to -M')
    %
end
%
disp(['running on ', int2str(poolsize), ' workers'])
%
wspace_filename = [wspace_filename, '_core', int2str(poolsize)];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of 2D_SS_Ising_pll_inner_v4_L', int2str(L), '_R1E', int2str(log10(REP))]);
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
    output = nan(q_max,12);
    %
else
    %
    output = nan(q_max, 10);
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
tStart = tic;
%
for q = 1:(q_max-1)
    %
    tic
    %
    neo_previous_all = cell(poolsize,1);
    hist_E_selected_S_vector_all = cell(poolsize,1);
    hist_WL_all = cell(poolsize,1);
    accept_counter_all = cell(poolsize,1);
    reject_counter_all = cell(poolsize,1);
    k_saved_all = cell(poolsize,1);
    k_all = cell(poolsize,1);
    sample_ratio_all = cell(poolsize,1);
    %
    parfor K_BIG = 1:poolsize % PARFOR
        %
        [neo_previous_all{K_BIG}, hist_E_selected_S_vector_all{K_BIG}, hist_WL_all{K_BIG}, accept_counter_all{K_BIG}, reject_counter_all{K_BIG}, k_saved_all{K_BIG}, k_all{K_BIG}, sample_ratio_all{K_BIG}] = function_flat_test_inner3_v53(q, N_atm, NN, skip, E_list, REP, poolsize, JDOS_aprox_frac, max_rw_steps, nnxpos, nnxneg, nnypos, nnyneg);
        %
    end
    %
    neo_previous = zeros(length(E_list), length(E_list)); % old new
    hits = JDOS_aprox_frac(:, q) > 0;
    %
    accept_counter = 0;
    reject_counter = 0;
    k_saved = 0;
    k = 0;
    sample_ratio = 0;
    %
    hist_WL = zeros(length(E_list),1);
    hist_E_selected_S_vector = zeros(length(E_list),1);
    %
    for K_BIG = 1:poolsize
        %
        neo_previous = neo_previous + neo_previous_all{K_BIG};
        hist_E_selected_S_vector = hist_E_selected_S_vector + hist_E_selected_S_vector_all{K_BIG};
        hist_WL = hist_WL + hist_WL_all{K_BIG};
        accept_counter = accept_counter + accept_counter_all{K_BIG};
        reject_counter = reject_counter + reject_counter_all{K_BIG};
        k_saved = k_saved + k_saved_all{K_BIG};
        k = k + k_all{K_BIG};
        sample_ratio = sample_ratio + sample_ratio_all{K_BIG};
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
    JDOS_aprox(:,q+1) = JDOS_aprox_frac(:,q+1) .* norm_factor(q+1);
    %
    q_timer = toc;
    %
    if L == 4
        %
        JDOS_erro(JDOS_exact_frac(:,q) ~= 0, q) = (JDOS_aprox_frac(JDOS_exact_frac(:,q) ~= 0,q) - JDOS_exact_frac(JDOS_exact_frac(:,q) ~= 0, q)) ./ JDOS_exact_frac(JDOS_exact_frac(:,q) ~= 0, q);
        JDOS_erro_total = JDOS_erro_total + sum(abs(JDOS_erro(:,q)));
        JDOS_erro_max(q,1) = max(abs(JDOS_erro(:,q)));
        %
        disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', num2str(nnz(hits)), ' | time per E pt: ', num2str(q_timer/nnz(hits)), ' secs | rw steps: ', int2str(k), ' | accept/reject: ', num2str(accept_counter/reject_counter), ' | avg. sample ratio: ', num2str(sample_ratio/poolsize), ' | JDOS_err_q: ', num2str(sum(abs(JDOS_erro(:,q)))), ' | cum JDOS_err: ', num2str(JDOS_erro_total)]);
        %
        output(q, 1:12) = [q, q_max, q_timer, nnz(hits), q_timer/nnz(hits), k, accept_counter/reject_counter, accept_counter, reject_counter, sample_ratio/poolsize, sum(abs(JDOS_erro(:,q))), JDOS_erro_total];
        %
    else
        %
        disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', num2str(nnz(hits)), ' | time per E pt: ', num2str(q_timer/nnz(hits)), ' secs | rw steps: ', int2str(k), ' | accept/reject: ', num2str(accept_counter/reject_counter), ' | avg. sample ratio: ', num2str(sample_ratio/poolsize)])
        %
        output(q, 1:10) = [q, q_max, q_timer, nnz(hits), q_timer/nnz(hits), k, accept_counter/reject_counter, accept_counter, reject_counter, sample_ratio/poolsize];
        %
    end
    %
end
%
q_timer = toc;
JDOS_aprox(:,q+1) = JDOS_aprox_frac(:,q+1) .* norm_factor(q+1);
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' | q: ', int2str(q+1), '/', int2str(q_max), ' time ', num2str(q_timer), ' secs']);
%
tEnd = toc(tStart);
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(tEnd), ' secs']);
%
if q_max == (L^2/2 + 1) || q_max == N_atm + 1
    %
    eval(['save ', wspace_filename, ' -v7.3'])
    %
end