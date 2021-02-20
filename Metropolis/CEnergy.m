function E = CEnergy(spins, J, L)
% Energy of a 2D Ising lattice with periodic boundary conditions.
E = 0;
aux = 1:L;
menos = circshift(aux, [0 -1]);
mais = circshift(aux, [0 +1]);
for i = 1:L
    for j = 1:L
        S = spins(i, j);
        WF = spins(mais(i), j) + spins(i, mais(j))...
            + spins(menos(i), j) + spins(i, menos(j));
        E = E - J * WF * S;
    end
end
E = E / 2;
end