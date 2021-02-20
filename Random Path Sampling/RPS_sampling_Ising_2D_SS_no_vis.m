clear
close all
%
L = 4; % linear size of system
REP = 1E6; % number of RPS sweeps // equivalent to (N_atm + 1)*REP spin configurations
%
NN = 4; % number of nearest neighbours
N_atm = L^2; % total number of spins
%
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = 1/2 * N_atm * NN : -4 : -1/2 * N_atm * NN; % possible energy values
%
% PRELIMINARY CALCULATION OF NEIGHBOUR TABLES
%
nnxpos = nan(L^2,1); % right neighbout in x
nnxneg = nan(L^2,1); % left neighbour in x
nnypos = nan(L^2,1); % up neighbour in y
nnyneg = nan(L^2,1); % down neighbour in z 
%
for i=1:L % loop through all row positions
    %
    for j=1:L % loop through all column positions
        %
        [nnxpos(j+(i-1)*L), nnxneg(j+(i-1)*L), nnypos(j+(i-1)*L), nnyneg(j+(i-1)*L)] = function_NN_list_2D_SS(L,i,j);
        %
    end
    %
end
%
E_all = nan(REP, length(M_list)); % full Energy matrix of all sweeps
E_all(:,1) = -1/2 * N_atm * NN; % energy of all spins pointing up
E_all(:,length(M_list)) = -1/2 * N_atm * NN; % energy of all spins pointing down
%
% RPS SWEEPS
%
RPS_timer = tic; % timer for RPS sampling
%
for k = 1:REP % loop through all requested RPS loops
    %
    S_vector = ones(N_atm, 1); % vector with spins
    SFV(:,1) = randperm(N_atm); % spin flip vector (sequence of spins to flip)
    %
    for q = 2:N_atm % loop through magnetization values
        %
        S_vector(SFV(q-1)) = -1; % flip the spin
        %
        E_new = - S_vector(SFV(q-1)) .* ( ...
            S_vector(nnxpos((SFV(q-1)))) + ...
            S_vector(nnxneg((SFV(q-1)))) + ...
            S_vector(nnypos((SFV(q-1)))) + ...
            S_vector(nnyneg((SFV(q-1))))); % energy of bonds to NN
        %
        E_all(k, q) = E_all(k, q-1) + 2*E_new; % build the energy matrix
        %
    end
    %
end
%
RPS_time = toc(RPS_timer); % register timer
disp(['RPS time ', num2str(RPS_time), ' seconds']); % display total RPS time
%
% energy histogram
%
hist_E = nan(numel(E_list), 2);
hist_E(:,1) = E_list;
%
for E_index = 1:numel(E_list)
    %
    hist_E(E_index,2) = sum(sum(E_all == E_list(E_index)));
    %
end
%
% magnetization histogram 
%
hist_M = nan(numel(M_list), 2);
hist_M(:,1) = M_list;
hist_M(:,2) = REP;
%
% (M,E) histogram
%
hist_timer = tic; % timer
%
hist_EM = nan(numel(E_list), numel(M_list));
%
for E_index = 1:numel(E_list)
    %
    for M_index = 1:numel(M_list)
        %
        hist_EM(E_index, M_index) = nnz(E_all(:, M_index) == E_list(E_index)) ;
        %
    end
    %
end
%
hist_time = toc(hist_timer); % register (E,M) histogram timer
disp(['(E,M) histogram time ', num2str(hist_time), ' seconds']); % display total histogram time
disp(['RPS + histogram time ', num2str(RPS_time + hist_time), ' seconds']); % display total time
%
% Histogram normalization to obtain the JDOS estimate
%
eval(['load norm_factor_Ising_Natm_',int2str(N_atm),'.mat'])
%
JDOS = nan(numel(E_list), numel(M_list)); 
%
for q = 1:N_atm
    %
    JDOS(:,q) = hist_EM(:,q)/REP * norm_factor(q);
    %
end