clear
close

% Metropolis Algorithm for the 3D Ising Model
% João Inácio, Sep. 17, 2020

    % Some variables
L = 4;              % Number of spins for each side of the lattice
% Use 16 as min
NSpins = L^3;       % Number of spins
NN = 6;             % Number of neighbours
J = 1;              % Interaction strenght between particles
T = 1.5:0.25:7.5;    % Temperature
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
    spins = ones(L, L, L);
    EConf = - 1/2  * J * NN * NSpins;
    
    % Magnetization and Energy matrix
    M = zeros(1, itMax);
    
    for it = 1:itMax
        i = randi([1 L]);
        j = randi([1 L]);
        k = randi([1 L]);
        S = spins(i, j, k);
        
        sumNei = spins(plus(i), j, k) + spins(minus(i), j, k)...
            + spins(i, plus(j), k) + spins(i, minus(j), k)...
            + spins(i, j, plus(k)) + spins(i, j, minus(k));
        dE = 2 * J * S * sumNei;
        ratio = exp(-dE * beta(nT));
        
        if dE <= 0 || rand() <= ratio
            spins(i, j, k) = - S;
        end
        
        M(it) = sum(sum(sum(spins)));
    end
    
    Mmean(nT) = mean(abs(M(end - (1E5 + 1):end)));
end
toc


figure(1)
plot(T, Mmean/NSpins, '-')
xlabel("Temperature")
ylabel("Mean Magnetization per Spin")












