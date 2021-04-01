clear
close all
%
% v72 - from Ising spin S 2D SS v7
%
% rng default
RandStream.setGlobalStream( RandStream.create('mt19937ar','seed',100) );
%
L = 4;
%
% q_max = 9;
%
dim = '2D';
lattice_type = 'SS';
neighbours = '1NN';
%
REP = 1E3; % number of desired configurations per (E,M) pair
skip = L^2; % 1 for no skip
%
% END OF USER INPUT
%
N_atm = L^2;
NN = 4;
%
if rem(length(-N_atm : 2 : N_atm),2) == 0
    %
    q_max = (length(-N_atm : 2 : N_atm))/2 - 1;
    %
else
    %
    q_max = (length(-N_atm : 2 : N_atm) + 1)/2 - 1; 
    %
end
%
disp(['Ising ', dim, ' ', lattice_type, ' ', neighbours])
disp(['L = ', int2str(L)])
disp(['REP = 1E', num2str(log10(REP))])
disp(['skip = ', int2str(skip)])
%
wspace_filename = ['workspace_nRPS_Ising_v72_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)];
JDOS_filename = ['JDOS_nRPS_Ising_v72_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', 'nRPS_Ising_v72_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)]);
%
eval(['load ./coefficients/coefficients_', int2str(N_atm), 'd', int2str(2),'.mat'])
eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])
%
output = nan(q_max, 5);
%
M_list(:,1) = -N_atm : 2 : N_atm;
E_list(:,1) = (- N_atm * NN /2) : 4 : (N_atm * NN /2); % possible energy values
%
JDOS_nRPS = zeros(length(E_list), length(M_list));
JDOS_nRPS(1,1) = 1;
JDOS_nRPS(1,length(M_list)) = 1;
%
t_total = tic;
%
% SCAN AT q = 1, ADD TO JDOS AT q = 2
[JDOS_nRPS] = ...
    function_Ising_scan_norm_correct_q1(N_atm, E_list, NN, NN_table, JDOS_nRPS);
%
% CALC JDOS at q = 2
JDOS_nRPS(:,2) = JDOS_nRPS(:,2) ./ sum(JDOS_nRPS(:,2)) .* norm_factor(2);
%
output(1,:) = [1, 0, 1, 0, 0];
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(1), '/', int2str(q_max)]);
%
% MAIN LOOP
%
for q = 2:q_max
    %
    q_timer = tic;
    %
    hist_nRPS = zeros(length(E_list), 1);
    hist_E_selected_nRPS = zeros(length(E_list), 1);
    %
    % RANDOM SPIN CONFIGURATION AT q
    [S_vector, E] = ...
        function_Ising_random_spin_config_at_q(N_atm, NN, q, NN_table);
    %
    hist_nRPS(E_list == E, 1) = hist_nRPS(E_list == E, 1) + 1;
    hist_E_selected_nRPS(E_list == E, 1) = hist_E_selected_nRPS(E_list == E, 1) + 1;
    %
    % SCAN
    [JDOS_nRPS] = ...
        function_Ising_scan_norm_correct(S_vector, E, E_list, NN_table, REP, JDOS_nRPS, q);
    %
    k = 1;
    %
    while min(hist_E_selected_nRPS(hist_E_selected_nRPS > 0)) < REP % CHECK FOR FULL CONFIG SET
        %
        % RANDOM WALK TRIAL STEP AT q
        [S_vector_new, E_new] = ...
            function_Ising_rw_step_at_q(S_vector, N_atm, E, NN_table);
        %
        % ACCEPT/REJECT WITH WL WEIGH FACTOR
        [S_vector, E, hist_nRPS] = ...
            function_Ising_WL_rw_criteria_correct(S_vector, E, S_vector_new, E_new, JDOS_nRPS(:,q), E_list, hist_nRPS);
        %
        % SCAN
        if hist_E_selected_nRPS(E_list == E, 1) < REP && rem(k,skip) == 0
            %
            [JDOS_nRPS] = ...
                function_Ising_scan_norm_correct(S_vector, E, E_list, NN_table, REP, JDOS_nRPS, q);            
            %
            hist_E_selected_nRPS(E_list == E, 1) = hist_E_selected_nRPS(E_list == E, 1) + 1;
            %
        end
        %
        k = k + 1;
        %
    end
    %
    JDOS_nRPS(:,q+1) = JDOS_nRPS(:,q+1) ./ sum(JDOS_nRPS(:,q+1)) .* norm_factor(q+1);
    %
    hits = nnz(JDOS_nRPS(:,q));
    q_timer = toc(q_timer);
    %
    output(q,:) = [q, q_timer, hits, q_timer/hits, k];
    disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ' int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', int2str(hits), ' | time per E pt: ', num2str(q_timer/hits), ' secs | rw steps: ', int2str(k)])
    %
end
%
t_total = toc(t_total);
%
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(t_total), ' seconds'])
eval(['save ', wspace_filename, '.mat -v7.3'])
eval(['save ', JDOS_filename, '.mat JDOS_nRPS E_list M_list -v7.3']) 
%
