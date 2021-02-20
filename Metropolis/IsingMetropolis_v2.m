clear
close

% Metropolis Algorithm for the Ising Model
% João Inácio, Sep. 4, 2020

    % Some variables
L = 50;            % Number of spins for each side of the lattice
% Use 16 as min
NSpins = L * L;     % Number of spins
nConfig = 2^NSpins; % Number of different spin configurations
J = 1;              % Interaction strenght between particles
T = 2.5;              % Temperature. We can change this value to observe some
% magnetic properties. If T > 2.26 the system is paramagnetic and for T <
% 2.26 the system is ferromagnetic
kB = 1;             % Boltzman constant in reduce units
beta = 1 / (T * kB);
itMax = 1E6;        % Maximum number of iterations for the method

% Vectors for the neighbours
aux = 1:L;
plus = circshift(aux,[0 -1]);
minus = circshift(aux,[0 +1]);

    % Metropolis Algorithm

% For this algorithm it's best to start with +1 spin matrix
spins = ones(L);
% spins = randSpins(L);

% Magnetization matrix
M = zeros(1, itMax);

tic
for it = 1:itMax
    i = randi([1 L]);
    j = randi([1 L]);
        
    sumNei = spins(plus(i), j) + spins(i, plus(j))...
            + spins(minus(i), j) + spins(i, minus(j));
    dE = 2 * J * spins(i, j) * sumNei;
    
    ratio = exp(-dE * beta);
    
    if dE <= 0 || rand() <= ratio
        spins(i, j) = - spins(i, j);
    end
    M(it) = sum(sum(spins));
end
toc

figure(1)
plot(1:itMax, M/NSpins)
title("Temperature = " + T + ", Lattice = " + L)
xlabel("Iterations")
ylabel("Magnetization")





