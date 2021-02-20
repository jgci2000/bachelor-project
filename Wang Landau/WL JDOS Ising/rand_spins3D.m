function spins = rand_spins3D(L)
% Generates a L*L*L matrix of random spins.

spins = zeros(L, L, L);

for i = 1:L
    for j = 1:L
        for k = 1:L
            spins(i, j, k) = sign(2 * rand() - 1);
        end
    end
end
end