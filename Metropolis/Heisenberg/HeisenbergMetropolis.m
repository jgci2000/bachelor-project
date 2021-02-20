clear
close

% Metropolis Algorithm for the 3D Ising Model
% João Inácio, Sep. 17, 2020

    % Some variables
L = 8;             % Number of spins for each side of the lattice
NSpins = L^3;       % Number of spins
NN = 6;             % Number of neighbours
J = 1;              % Interaction strenght between particles
T = 0:0.5:5;    % Temperature. We can change this value to observe some
% magnetic properties. If T > 2.26 the system is paramagnetic and for T <
% 2.26 the system is ferromagnetic
kB = 1;             % Boltzman constant in reduce units
beta = ones(1, length(T)) ./ (T * kB);
itMax = 1E6;        % Maximum number of iterations for the method

% Vectors for the neighbours
aux = 1:L;
plus = circshift(aux, [0 -1]);
minus = circshift(aux, [0 +1]);

    % Metropolis Algorithm

Mmean = zeros(1, length(T));

tic
for nT = 1:length(T)
    % For this algorithm it's best to start with +1 spin matrix
    % spins{1} = spinsX for each spin in an L*L*L lattice
    spins = cell(1, 3);
    spins{1} = zeros(L, L, L);
    spins{2} = zeros(L, L, L);
    spins{3} = ones(L, L, L);
    
    % Magnetization and Energy matrix
    M = zeros(1, itMax);
    
    for it = 1:itMax
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
        ratio = exp(-dE * beta(nT));
        
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
    
    Mmean(nT) = mean(abs(M(end - (1E5 + 1):end)));
end
toc


figure(1)
plot(T, Mmean/NSpins, '-')
xlabel("Temperature")
ylabel("Mean Magnetization per Spin")












