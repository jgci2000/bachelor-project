function M = Heisenberg(L, T, it_max)

    % Some variables
NSpins = L^3;       % Number of spins
NN = 6;             % Number of neighbours
J = 1;              % Interaction strenght between particles
kB = 1;             % Boltzman constant in reduce units
beta = 1 / (T * kB);

% Vectors for the neighbours
aux = 1:L;
plus = circshift(aux, [0 -1]);
minus = circshift(aux, [0 +1]);

    % Metropolis Algorithm

% For this algorithm it's best to start with +1 spin matrix
% spins{1} = spinsX for each spin in an L*L*L lattice
spins = cell(1, 3);
spins{1} = zeros(L, L, L);
spins{2} = zeros(L, L, L);
spins{3} = ones(L, L, L);

% Magnetization and Energy matrix
M = zeros(1, it_max);
    
for it = 1:it_max
    i = randi([1 L]);
    j = randi([1 L]);
    k = randi([1 L]);

    Sx = spins{1}(i, j, k);
    Sy = spins{2}(i, j, k);
    Sz = spins{3}(i, j, k);

    % Marsaglia method (new spin components)
    x1 = 1;
    x2 = 1;
    while (x1 ^ 2 + x2 ^ 2 >= 1)
        x1 = - 1 + 2 * rand();
        x2 = - 1 + 2 * rand();
    end

    NSx = 2 * x1 * sqrt(1 - x1 ^ 2 - x2 ^ 2);
    NSy = 2 * x2 * sqrt(1 - x1 ^ 2 - x2 ^ 2);
    NSz = 1 - 2 * (x1 ^ 2 + x2 ^ 2);

    sumNeix = spins{1}(plus(i), j, k) + spins{1}(minus(i), j, k)...
        + spins{1}(i, plus(j), k) + spins{1}(i, minus(j), k)...
        + spins{1}(i, j, plus(k)) + spins{1}(i, j, minus(k));
    sumNeiy = spins{2}(plus(i), j, k) + spins{2}(minus(i), j, k)...
        + spins{2}(i, plus(j), k) + spins{2}(i, minus(j), k)...
        + spins{2}(i, j, plus(k)) + spins{2}(i, j, minus(k));
    sumNeiz = spins{3}(plus(i), j, k) + spins{3}(minus(i), j, k)...
        + spins{3}(i, plus(j), k) + spins{3}(i, minus(j), k)...
        + spins{3}(i, j, plus(k)) + spins{3}(i, j, minus(k));
    dE = - J * ((NSx - Sx) * sumNeix + (NSy - Sy) * sumNeiy + (NSz - Sz) * sumNeiz);
    ratio = exp(-dE * beta);

    if dE <= 0 || rand() <= ratio
        spins{1}(i, j, k) = NSx;
        spins{2}(i, j, k) = NSy;
        spins{3}(i, j, k) = NSz;
    end

    Mx = sum(sum(sum(spins{1})));
    My = sum(sum(sum(spins{2})));
    Mz = sum(sum(sum(spins{3})));

    M(it) = sqrt(Mx ^ 2 + My ^ 2 + Mz ^ 2);
end

end

