clear
close all
%
% v1 - first implementation, works for Npos = 2
% v2 - attempt to generalize for Npos > 2 % OK, need to fix JDOS calcs
% v3 - fixing JDOS calcs % works for Npos = 3,4 + L = 4,8, but not L = 2?
% v4 - fix for L = 2? % no, for some reason, worse for L = 4! real fix = hard
% v5 - with Npos_sum_config tables, try correct calculations | WORKS (for Npos=3)
% v6 - trying to fix for Npos > 3 (yes!)
% v7 - cleanup and QoL (output)
%
rng default
%
L = 4;
Npos = 4; % (2*S)+1
%
% q_max = 11; %ceil((length(-L^2 : 2/(Npos-1) : L^2/2) + 1)/2);
%
dim = '2D';
lattice_type = 'SS';
neighbours = '1NN';
%
REP = 1E3; % number of desired configurations per (E,M) pair
skip = 1; %L^2*(Npos-1); % 1 for no skip
%
% END OF USER INPUT
%
N_atm = L^2;
NN = 4;
%
if rem(length(-double(N_atm) : 2/double(Npos-1) : double(N_atm)),2) == 0
    %
    q_max = (length(-double(N_atm) : 2/double(Npos-1) : double(N_atm)))/2 - 1; %length(-N_atm : 2/(Npos-1) : N_atm);
    %
else
    %
    q_max = (length(-double(N_atm) : 2/double(Npos-1) : double(N_atm)) + 1)/2 - 1; %length(-N_atm : 2/(Npos-1) : N_atm);
    %
end
%
disp(['Spin S Ising ', dim, ' ', lattice_type, ' ', neighbours])
disp(['Npos = ', int2str(Npos)])
disp(['L = ', int2str(L)])
disp(['REP = 1E', num2str(log10(REP))])
disp(['skip == ', int2str(skip)])
%
wspace_filename = ['workspace_nRPS_spinS_Ising_v7_Npos', int2str(Npos), '_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)];
JDOS_filename = ['JDOS_nRPS_spinS_Ising_v7_Npos', int2str(Npos), '_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', 'nRPS_spinS_Ising_v7_Npos', int2str(Npos), '_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)]);
%
eval(['load ./coefficients/coefficients_', int2str(N_atm), 'd', int2str(Npos),'.mat'])
eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])
eval(['load ./Npos_sum_configs/Npos_sum_configs_Npos', int2str(Npos), '_N_atm', int2str(N_atm), '.mat'])
%
output = nan(q_max, 5);
%
Z_spin_values(:,1) = -(Npos-1) : 2 : (Npos-1); % always integers
%
M_list(:,1) = -N_atm*max(Z_spin_values) : 2 : N_atm*max(Z_spin_values);
E_list(:,1) = max(Z_spin_values(:,1)).^2*(- N_atm * NN ./2) : 4 : max(Z_spin_values(:,1)).^2*(N_atm * NN ./2); % possible energy values
%
JDOS_nRPS = zeros(length(E_list), length(M_list));
JDOS_nRPS(1,1) = 1;
JDOS_nRPS(1,length(M_list)) = 1;
%
% CREATE JDOS_nRPS_correct
JDOS_nRPS_correct = cell(length(M_list),1);
%
for q = 1:length(M_list)
    %
    JDOS_nRPS_correct{q,1} = zeros(length(E_list), length(Npos_sum_configs{q,1}(:,1)));
    %
end
%
JDOS_nRPS_correct{1,1}(1,1) = 1;
JDOS_nRPS_correct{length(M_list),1}(1,1) = 1;
%
t_total = tic;
%
% SCAN AT q = 1, ADD TO JDOS_correct q = 2+
[JDOS_nRPS_correct] = ...
    function_scan_norm_correct_q1(Npos, N_atm, E_list, NN, NN_table, Z_spin_values, Npos_sum_configs, JDOS_nRPS_correct);
%
% CALC JDOS at q = 2
JDOS_nRPS(:,2) = sum(JDOS_nRPS_correct{2,1},2) ./ sum(sum(JDOS_nRPS_correct{2,1})) .* norm_factor(2);
%
output(1,:) = [1, 0, 1, 0, 0];
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(1), '/', int2str(q_max)]);
%
% MAIN LOOP
for q = 2:q_max
    %
    q_timer = tic;
    %
    hist_nRPS = zeros(length(E_list), length(Npos_sum_configs{q,1}(:,1)));
    hist_E_selected_nRPS = zeros(length(E_list), length(Npos_sum_configs{q,1}(:,1)));
    %
    % RANDOM SPIN CONFIGURATION AT q
    [S_vector, SPM, E, index_Npos_sum_vector] = ...
        function_random_spin_config_at_q(N_atm, NN, Z_spin_values, q, Npos, NN_table, Npos_sum_configs{q,1});
    %
    hist_nRPS(E_list == E, index_Npos_sum_vector) = hist_nRPS(E_list == E, index_Npos_sum_vector) + 1;
    hist_E_selected_nRPS(E_list == E, index_Npos_sum_vector) = hist_E_selected_nRPS(E_list == E, index_Npos_sum_vector) + 1;
    %
    % SCAN
    [JDOS_nRPS_correct] = ...
        function_scan_norm_correct(S_vector, SPM, E, Npos, E_list, NN_table, Z_spin_values, REP, Npos_sum_configs, JDOS_nRPS_correct, q);
    %
    k = 1;
    %
    while min(hist_E_selected_nRPS(hist_E_selected_nRPS > 0)) < REP % CHECK FOR FULL CONFIG SET
        %
        % RANDOM WALK TRIAL STEP AT q
        [S_vector_new, SPM_new, E_new, index_Npos_sum_vector_new] = ...
            function_rw_step_at_q(S_vector, SPM, N_atm, Npos, E, NN_table, Z_spin_values, Npos_sum_configs{q,1});
        %
        % ACCEPT/REJECT WITH WL WEIGH FACTOR
        [S_vector, SPM, E, hist_nRPS, index_Npos_sum_vector] = ...
            function_WL_rw_criteria_correct(S_vector, SPM, E, S_vector_new, SPM_new, E_new, JDOS_nRPS_correct{q,1}, E_list, hist_nRPS, index_Npos_sum_vector, index_Npos_sum_vector_new);
        %
        % SCAN
        if hist_E_selected_nRPS(E_list == E, index_Npos_sum_vector) < REP && rem(k,skip) == 0
            %
            [JDOS_nRPS_correct] = ...
                function_scan_norm_correct(S_vector, SPM, E, Npos, E_list, NN_table, Z_spin_values, REP, Npos_sum_configs, JDOS_nRPS_correct, q);
            %
            hist_E_selected_nRPS(E_list == E, index_Npos_sum_vector) = hist_E_selected_nRPS(E_list == E, index_Npos_sum_vector) + 1;
            %
        end
        %
        k = k + 1;
        %
    end
    %
    JDOS_nRPS_correct{q+1,1} = JDOS_nRPS_correct{q+1,1} ./ sum(sum(JDOS_nRPS_correct{q+1,1})) .* norm_factor(q+1);
    JDOS_nRPS(:,q+1) = sum(JDOS_nRPS_correct{q+1,1},2) ./ sum(sum(JDOS_nRPS_correct{q+1,1})) .* norm_factor(q+1);
    %
    hits = nnz(JDOS_nRPS_correct{q,1});
    q_timer = toc(q_timer);
    %
    output(q,:) = [q, q_timer, hits, q_timer/hits, k];
    disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', int2str(hits), ' | time per E pt: ', num2str(q_timer/hits), ' secs | rw steps: ', int2str(k)])
    %
end
%
t_total = toc(t_total);
%
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(t_total), ' seconds'])
eval(['save ', wspace_filename, '.mat -v7.3'])
eval(['save ', JDOS_filename, '.mat JDOS_nRPS E_list M_list -v7.3']) 
%