clear
close all
%
% WL estimate of DOS at a given M (from q value)
% norm_factor in log, for calculations of large systems >32x32
% main outputs:   
%   E_list - list of energy values
%   WL_DOS_f_norm - normalized DOS
%   WL_log_DOS_f_norm - log of normalized DOS
%   output - (f value, time in seconds, number of steps)
%
rng default
%
L = 4;
q = L^2/2 + 1;
f = exp(1);
p = 0.95;
f_min = 1 + 1E-6;
%
check_flat = 1E4;
max_rw_steps = 1E14;
%
% END OF USER INPUT
%
wspace_filename = ['workspace_WL_sqrt_f_log_norm_factor_v1_L', int2str(L), '_q', int2str(q), '.mat'];
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', wspace_filename]);
%
NN = 4;
N_atm = L^2;
E_list(:,1) = -1 * N_atm * NN /2 : +4 : N_atm * NN /2; % possible energy values
%
eval(['load norm_factor_Ising_Natm_', int2str(N_atm), '_log.mat']);
%
% LATTICE AND NEIGHBOUR TABLES
[nnxpos, nnxneg, nnypos, nnyneg] = function_neighbours_2D_SS(L);
%
tic
tstart = tic;
%
WL_log_DOS_dynamic = zeros(length(E_list), 1);
%
hist_WL = zeros(length(E_list),1);
hist_E_selected_S_vector = zeros(length(E_list),1);
%
% RANDOM CONFIGURATION AT q
[S_vector_WL, E_WL_old] = function_random_config_at_q_2D_SS(N_atm, NN, q, nnxpos, nnxneg, nnypos, nnyneg);
%
hist_WL(E_list == E_WL_old) = hist_WL(E_list == E_WL_old) + 1;
WL_log_DOS_dynamic(E_list == E_WL_old, 1) = WL_log_DOS_dynamic(E_list == E_WL_old, 1) + log(f);
%
f_temp = f;
%
while f_temp > f_min
    %
    f_temp(length(f_temp)+1,1) = sqrt(f_temp(length(f_temp)));
    %
end
%
WL_log_DOS_f = zeros(length(E_list), length(f_temp));
%
output = nan(length(f_temp),3);
output(:,1) = f_temp(:,1)-1;
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | f = 1 + ', num2str(f-1, '%.4e'), ' start']);
%
figure(1)
subplot(1,2,1)
plot(E_list, hist_WL, 'ko-')
subplot(1,2,2)
plot(E_list, WL_log_DOS_dynamic, 'ko-')
%
t_start = tic;
t_f_start = tic;
%
k_saved = 0;
f_counter = 1;
%
for k = 2:max_rw_steps
    %
    if rem(k, check_flat) == 0
        %
        subplot(1,2,1)
        plot(E_list, hist_WL, 'ko-')
        subplot(1,2,2)
        plot(E_list, WL_log_DOS_dynamic, 'ko-')
        drawnow
        %
        if min(hist_WL(hist_WL > 0)) >= p * mean(hist_WL(hist_WL > 0)) && max(hist_WL(hist_WL > 0)) <= (1 + (1-p)) * mean(hist_WL(hist_WL > 0))
            %
            t_f_end = toc(t_f_start);
            disp([int2str(f_counter), '/', int2str(length(f_temp)), ' | ', datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | f = 1 + ', num2str(f-1, '%.4e'), ' end | time: ', num2str(t_f_end), ' | rw steps: ', int2str(k - k_saved)]);
            output(f_counter,2) = t_f_end;
            output(f_counter,3) = k - k_saved;
            k_saved = k;
            WL_log_DOS_f(:, f_counter) = WL_log_DOS_dynamic;
            %
            if f < f_min
                %
                break
                %
            end
            %
            hist_WL = zeros(length(E_list),1);
            f = sqrt(f);
            f_counter = f_counter + 1;
            %
            disp([int2str(f_counter), '/', int2str(length(f_temp)), ' | ', datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | f = 1 + ', num2str(f-1, '%.4e'), ' start']);
            t_f_start = tic;
            %
        end
        %
    end
    %
    [S_vector_WL, E_WL_old, WL_log_DOS_dynamic, hist_WL] = ...
        function_WL_RW_step_no_scan(S_vector_WL, E_WL_old, nnxpos, nnxneg, nnypos, nnyneg, E_list, hist_WL, WL_log_DOS_dynamic, f);
    %
end
%
t_total = toc(t_start);
%
WL_log_DOS_f_norm = zeros(length(E_list), length(f_temp));
WL_DOS_f_norm = zeros(length(E_list), length(f_temp));
%
for f_counter = 1:length(f_temp)
    %
    hits = WL_log_DOS_f(:, f_counter) > 0;
    %
    list_hits = find(hits);
    %
    sum_log_DOS_temp = WL_log_DOS_f(list_hits(1), f_counter);
    %
    for index_hits = 2:length(list_hits)
        %
        sum_log_DOS_temp = sum_log_DOS_temp + log(1 + exp(WL_log_DOS_f(list_hits(index_hits),  f_counter) - sum_log_DOS_temp));
        %
    end
    %
    WL_log_DOS_f_norm(hits, f_counter) = WL_log_DOS_f(hits, f_counter) - sum_log_DOS_temp + norm_factor_log(q);
    WL_DOS_f_norm(WL_log_DOS_f_norm > 0) = exp(WL_log_DOS_f_norm(WL_log_DOS_f_norm > 0));
    %
end
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | total runtime: ', num2str(t_total), ' secs.'])
%
eval(['save ', wspace_filename, ' -v7.3']) 
