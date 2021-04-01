function E = comp_energy3D(spins, J, L)
% Energy of a 2D Ising lattice with periodic boundary conditions.

E = 0;

aux = 1:L;
menos = circshift(aux, [0 -1]);
mais = circshift(aux, [0 +1]);

for i = 1:L
    for j = 1:L
        for k = 1:L
            s = spins(i, j, k);
            sum_nei = spins(mais(i), j, k) + spins(i, mais(j), k)...
                + spins(menos(i), j, k) + spins(i, menos(j), k)...
                + spins(i, j, mais(k)) + spins(i, j, menos(k));
            
            E = E - J * sum_nei * s;
        end
    end
end

E = E / 2;
end