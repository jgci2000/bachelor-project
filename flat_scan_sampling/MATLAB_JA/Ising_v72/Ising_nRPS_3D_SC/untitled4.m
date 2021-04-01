RandStream.setGlobalStream( RandStream.create('mt19937ar','seed',100) );
% rand(100,1)

L = 4;

dim = '3D';
lattice_type = 'SC';
neighbours = '1NN';

REP = 1E3; % number of desired configurations per (E,M) pair
skip = L^2; % 1 for no skip

N_SPINS = L^3;
NN = 4;

eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])

M_list(:, 1) = -N_SPINS : 2 : N_SPINS;
E_list(:, 1) = (- N_SPINS * NN /2) : 4 : (N_SPINS * NN /2); % possible energy values

spins = zeros(1, N_SPINS);

for i = 1:N_SPINS
    spins(i) = 1;
    if rand() < 0.75
        spins(i) = -1;
    end
end

neo_previous = zeros(length(E_list), length(E_list)); % old new

E = 0;
for i = 1:N_SPINS
    sum_nei = 0;
    for a = 1:NN
        sum_nei = sum_nei + spins(NN_table(i, a));
    end
    
    E = E - spins(i) * sum_nei;
end
E_WL_old = E / 2;

pos_scan = find(spins == 1);
E_old = E_WL_old;
delta_E = spins .* (spins(NN_table(:, 1)) + spins(NN_table(:, 2)) + spins(NN_table(:, 3)) + ...
    + spins(NN_table(:, 4)) + spins(NN_table(:, 5)) + spins(NN_table(:, 6)));
E_new = E_old + 2*delta_E(pos_scan);
for i = 1:length(E_new)
    neo_previous(E_list == E_old, :) = neo_previous(E_list == E_old, :) + sum(E_list == E_new(i),2)';
end
