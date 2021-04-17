function M = comp_magnetization3D(spins, L)
% Magnetization of a 2D Ising lattice.

M = 0;

for i = 1:L
    for j = 1:L
        for k = 1:L
            M = M + spins(i, j, k);
        end
    end
end

end

