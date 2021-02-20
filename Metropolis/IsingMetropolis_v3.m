clear
close

% Metropolis Algorithm for the Ising Model
% João Inácio, Sep. 4, 2020

    % Some variables
L = 4;            % Number of spins for each side of the lattice
% Use 16 as min
NSpins = L^2;       % Number of spins
J = 1;              % Interaction strenght between particles
T = 1.5:0.1:3.5;    % Temperature. We can change this value to observe some
% magnetic properties. If T > 2.26 the system is paramagnetic and for T <
% 2.26 the system is ferromagnetic
kB = 1;             % Boltzman constant in reduce units
beta = ones(1, length(T)) ./ (T * kB);
itMax = 1E6;        % Maximum number of iterations for the method

% Vectors for the neighbours
aux = 1:L;
plus = circshift(aux,[0 -1]);
minus = circshift(aux,[0 +1]);

    % Metropolis Algorithm

Mmean = zeros(1, length(T));
Emean = zeros(1, length(T));

tic
for k = 1:length(T)
    % For this algorithm it's best to start with +1 spin matrix
    spins = ones(L);
    EConf = - 2 * NSpins;
    
    % Magnetization and Energy matrix
    M = zeros(1, itMax);
    E = zeros(1, itMax);
    
    for it = 1:itMax
        i = randi([1 L]);
        j = randi([1 L]);
        
        sumNei = spins(plus(i), j) + spins(i, plus(j))...
            + spins(minus(i), j) + spins(i, minus(j));
        dE = 2 * J * spins(i, j) * sumNei;
        ratio = exp(-dE * beta(k));
        
        if dE <= 0 || rand() <= ratio
            spins(i, j) = - spins(i, j);
            EConf = dE + EConf;
        end
        
        M(it) = sum(sum(spins));
        E(it) = EConf;
    end
    
    Mmean(k) = mean(abs(M(end - (1E5 + 1):end)));
    Emean(k) = mean(E(end - (1E5 + 1):end));
end
toc

% Compute F
% F = - kB * T * ln(Z) and Z = sum(exp(- E  / beta))
Z = exp(- Emean ./ beta);
F = - kB * T .* log(Z);

figure(1)
subplot(2, 1, 1)
plot(T, Mmean/NSpins, '-')
xlabel("Temperature")
ylabel("Mean Magnetization per Spin")

subplot(2, 1, 2)
plot(Mmean/NSpins, F/NSpins, 'o')
xlabel("Free Energy")
ylabel("Mean Magnetization")










