clear
close all
%
% v1 - first version
%
% rng default
%
L = 8;
q_max = L^2/2 + 1;
%
max_rw_steps = 1E12;
skip = L^2; % 0 for no skip
REP = 1E5;
slice = 1E3;
%
NN = 4;
N_atm = L^2;
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = -1/2 * N_atm * NN : +4 : 1/2 * N_atm * NN; % possible energy values
%
wspace_filename = ['workspace_ft_v48_parallel_L', int2str(L), 'R1E', int2str(log10(slice)),'x',  '1E', int2str(log10(floor(REP/slice))), '.mat'];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of 2D_SS_Ising_v48_L' int2str(L), 'R1E', int2str(log10(slice)),'x',  '1E', int2str(log10(floor(REP/slice)))]);
disp(['output workspace filename: ', wspace_filename])
%
disp(['L = ', int2str(L)])
disp(['REP = 1E', num2str(log10(REP))])
disp(['slice = 1E', num2str(log10(slice))])
disp(['skip = ', int2str(skip)])
%
eval(['load norm_factor_Ising_Natm_', int2str(N_atm)])
%
JDOS_aprox_frac = cell(floor(REP/slice),1);
output = cell(floor(REP/slice),1);
%
tic
t_start = tic;
%
parfor BIG_K = 1:floor(REP/slice)
    %
    [JDOS_aprox_frac{BIG_K}, output{BIG_K}] = function_flat_test_v48(L, q_max, skip, max_rw_steps, slice, BIG_K);
    %
end
%
toc
t_end = toc(t_start);
%
JDOS_aprox = cell(floor(REP/slice),1);
avg_JDOS_aprox = zeros(length(E_list), length(M_list));
%
for BIG_K = 1:floor(REP/slice)
    %
    for q = 1:(N_atm+1)
        %
        JDOS_aprox{BIG_K}(:,q) = JDOS_aprox_frac{BIG_K}(:,q) .* norm_factor(q);
        %
    end
    %
    avg_JDOS_aprox = avg_JDOS_aprox + JDOS_aprox{BIG_K};
    %
end
%
avg_JDOS_aprox = avg_JDOS_aprox ./ floor(REP/slice);
%
DOS_slice = zeros(floor(REP/slice),1);
DOS_zerozero = zeros(floor(REP/slice),1);
DOS_chkbrd = zeros(floor(REP/slice),1);
%
for BIG_K = 1:floor(REP/slice)
    %
    DOS_slice(BIG_K,1) = JDOS_aprox{BIG_K}(5,9);
    DOS_zerozero(BIG_K,1) = JDOS_aprox{BIG_K}(9,9);
    DOS_chkbrd(BIG_K,1) = JDOS_aprox{BIG_K}(17,9);
    %
end
%
eval(['save ', wspace_filename, ' -v7.3'])
%