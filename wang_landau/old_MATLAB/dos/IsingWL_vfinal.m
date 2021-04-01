clear
close
clc

% Wang Landau sampling for 2D Ising Model
% João Inácio, Aug. 19, 2020

    % Constants
L = 2;              % Number of spins for each side of the lattice.
NSpins = L * L;     % Number of spins.
NConfig = 2^NSpins; % Number of different spin configurations
J = 1;              % Interaction strenght between particles.

    % Neighbours Vectors
aux = 1:L;
menos = circshift(aux, [0 -1]);
mais = circshift(aux, [0 +1]);

    % Wang Landau Sampling
f = exp(1);     % Modification factor.
p = 0.95;       % Flatness -> min(hist) > avg(hist)*p.

% All of the possible energies. EMin + 4 and EMax -4 are not possible.
temp = 4 * (0:NSpins) - 2 * NSpins;
energies(1) = temp(1);
energies(2:length(temp) - 3) = temp(3:end-2);
energies(length(temp) - 2) = temp(end);
NE = length(energies);

% Random configuration
spins = randSpins(L);
E = CEnergy(spins, J, L);
idxE = find(energies == E);

% Desnsity of States Vector
g = ones(1, NE);
% Working with the ln(g(E)) is better because the values are too high.
lngE = log(g);

% At the end of 10000 MCSweeps check if the histogram is flat.
MCSweeps = 0;   % Reset the sweeps counter.

% histogram
hist = zeros(1, NE);
tic
while f > 1 + 1E-8
    % Each for loop is 1 MCSweep
    for ni = 1:NSpins
        % Select a random spin
        i = randi([1 L]);
        j = randi([1 L]);
        S = spins(i, j);
        
        % dE = 2 * J * S * sum(neighbours) and dE = ENew - E, so
        sumNei = spins(mais(i), j) + spins(i, mais(j))...
            + spins(menos(i), j) + spins(i, menos(j));
        ENew = E + 2 * J * S * sumNei;
        idxEN = find(energies == ENew);
        
        % Flip the spin
        ratio = exp(lngE(idxE) - lngE(idxEN));
        P = min(ratio, 1);
        
        if P > rand()
            spins(i, j) = -S;
            E = ENew;
            idxE = idxEN;
        end
        
        % Update the histogram and g(E)
        hist(idxE) = hist(idxE) + 1;
        lngE(idxE) = lngE(idxE) + log(f);
    end
    
    MCSweeps = MCSweeps + 1;
    
    % Check the flatness of the histogram each 10000 setps.
    if mod(MCSweeps, 10000) == 0
        avgH = mean(hist);
        minH = min(hist);
        hist
        if minH > avgH * p
            hist = zeros(1, NE);
            fprintf("%d: the histogram is flat. Min: %.0f Avg: %.2f f: %.8f\n", MCSweeps, minH, avgH, f);
            f = f^(1/2);
            % Reset the counter.
            MCSweeps = 0;
        end
    end
end
toc

% Normalize the DOS
% lngE is the normalized Density of States.
lngE = lngE - lngE(energies == -2 * NSpins) + log(2);

% Get the actual DOS
g = exp(lngE);

figure(2)
plot(energies/NSpins, g)




