clear
close all
%
% v1 - first implementation
% v2 - generalization for NN (no change in output)
% v3 - generalization for E and WL functions (no change in output)
%
rng default
%
L = 8;
Npos = 2; % (2*S)+1
%
dim = '2D';
lattice_type = 'SS';
neighbours = '1NN';
%
REP = 1E5; % when to check for histogram flatness
max_rw_steps = 1E14; % maximum rw steps 
p = 0.95; % flatness criteria
log_f_start = 1E0; % starting log value of f
log_f_min = 1E-6; % target low value of f
%
% END OF USER INPUT
%
N_atm = L^2;
%
wspace_filename = ['workspace_WL_spinS_Ising_dev_Npos', int2str(Npos), '_2D_SS_L', int2str(L), '_p_', num2str(p), '_log_f_start_1E', int2str(log10(log_f_start)), '_log_f_min_1E', int2str(log10(log_f_min)), '_REP_1E', int2str(log10(REP))];
JDOS_filename = ['JDOS_WL_spinS_Ising_dev_Npos', int2str(Npos), '_2D_SS_L', int2str(L), '_p_', num2str(p), '_log_f_start_1E', int2str(log10(log_f_start)), '_log_f_min_1E', int2str(log10(log_f_min)), '_REP_1E', int2str(log10(REP))];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', wspace_filename]);
%
disp(['Npos = ', int2str(Npos)])
disp(['REP = 1E', num2str(log10(REP))])
disp(['log f_start = 1E', num2str(log10(log_f_start))])
disp(['log f_min = 1E', num2str(log10(log_f_min))])
disp(['p = ', num2str(p)])
%
eval(['load ./coefficients/coefficients_', int2str(N_atm), 'd', int2str(Npos),'.mat'])
eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])
%
Z_spin_values(:,1) = (Npos-1) : -2 : -(Npos-1); % always integers
%
NN = 4;
M_list(:,1) = -N_atm*max(Z_spin_values) : 2 : N_atm*max(Z_spin_values);
E_list(:,1) = max(Z_spin_values(:,1)).^2*(- N_atm * NN ./2) : 4 : max(Z_spin_values(:,1)).^2*(N_atm * NN ./2); % possible energy values
%
% PREALLOCATE JDOS
log_JDOS_WL = zeros(length(E_list), length(M_list));
%
log_f = log_f_start;
%
t_start = tic;
%
while log_f > (1/2*log_f_min)
    %
    % ZERO HISTOGRAM
    hist_WL = zeros(length(E_list), length(M_list));
    %
    % RANDOM CONFIGURATION
    SPM_vector = randi(Npos,[N_atm 1]);
    S_vector = Z_spin_values(SPM_vector);
    %
    % CALC E AND M_z
    E = function_Energy_Ising_2D_SS_neo_dev(L, S_vector, NN_table);
    M_z = sum(S_vector(:,1));
    %
    % UPDATE HISTOGRAM
    hist_WL(E_list == E, M_list == M_z) = 1;    
    %
    for k = 2:max_rw_steps
        %
        if rem(k, REP) == 0
            %
            % CHECK HISTOGRAM FLATNESS
            if min(hist_WL(hist_WL > 0)) >= p * mean(hist_WL(hist_WL > 0)) && max(hist_WL(hist_WL > 0)) <= (1 + (1-p)) * mean(hist_WL(hist_WL > 0))
                %
                disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | log f = ', num2str(log_f)])
                hist_WL = zeros(length(E_list), length(M_list));
                log_f = 1/2 * log_f;
                %
                break
                %
            end
            %
        end
        %
        % RANDOM STEP, WL CRITERIA, CALC E AND M_z, UPDATE HISTOGRAM AND log_JDOS
        [E, M_z, S_vector, hist_WL, log_JDOS_WL] = ...
            function_rw_step_WL_spinS_Ising_dev(E, M_z, S_vector, hist_WL, N_atm, Npos, Z_spin_values, log_JDOS_WL, log_f, E_list, M_list, NN_table);
        %
    end
    %
end
%
% CALC NORMALIZED JDOS
%
JDOS_WL = zeros(length(E_list), length(M_list));
%
for q = 1:length(M_list)
    %
    hits = log_JDOS_WL(:,q) > 0;
    %
    JDOS_WL(hits,q) = log_JDOS_WL(hits, q) - max(log_JDOS_WL(hits, q));
    JDOS_WL(hits, q) = exp(JDOS_WL(hits, q)) ./ sum(exp(JDOS_WL(hits, q))) * norm_factor(q);
    %
end
%
t_end = toc(t_start);
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | finished, time = ', num2str(t_end)])
%
eval(['save ', wspace_filename, '.mat -v7.3'])
eval(['save ', JDOS_filename, '.mat JDOS_WL E_list M_list -v7.3'])