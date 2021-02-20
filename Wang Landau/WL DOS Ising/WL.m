clear
close
clc


J = 1;
L = 2;
NSpins = L * L;
flatness = 0.9;

% Neighbours Vectors
aux = 1:L;
menos = circshift(aux, [0 -1]);
mais = circshift(aux, [0 +1]);

% Wang Landau algorithm
% All of the possible energies. EMin + 4 and EMax -4 are not possible.
temp = 4 * (0:NSpins) - 2 * NSpins;
Energies(1) = temp(1);
Energies(2:length(temp) - 3) = temp(3:end-2);
Energies(length(temp) - 2) = temp(end);

% Random configuration
spins = randS(L);
E = CEnergy(spins, J, L);
idxE = find(Energies == E);


% Natural logarithm of the density of states
lngE = zeros(1, length(Energies));
% Histogram
Hist = zeros(1, length(Energies));

% Modification factor
lnf = 1;

MCSweep = 0;

tic
while exp(lnf) > 1 + 1E-8
    for it = 1:NSpins
        i = randi([1 L]);
        j = randi([1 L]);
        
        % Compute the new energy
        S = spins(i, j);
        WF = spins(mais(i), j) + spins(i, mais(j))...
            + spins(menos(i), j) + spins(i, menos(j));
        ENew = E + 2 * J * S * WF;  % dE = 2 * J * S * WF, so...
        idxEN = find(Energies == ENew);
        
        % Flip the spin
        P = exp(lngE(idxE) - lngE(idxEN));
        if P > rand()
            spins(i, j) = -S;
            E = ENew;
            idxE = idxEN;
        end
        
        % Update the histogram and lngE
        Hist(idxE) = Hist(idxE) + 1;
        lngE(idxE) = lngE(idxE) + lnf;
    end
    
    MCSweep = MCSweep + 1;
    
    % Check the flatness each 10000 MCSweeps
    if MCSweep == 10000
        MCSweep = 0;
        
        avgH = mean(Hist);
        minH = min(Hist);
        
        if minH > avgH * flatness
            Hist = zeros(1, length(Hist));
            lnf = lnf / 2;
            fprintf("%d: the histogram is flat. Min: %f Avg: %f f: %f\n", MCSweep, minH, avgH, exp(lnf));
        end
    end
end
toc

% Normalize the DOS
if lngE(end) < lngE(1)
    lgC = lngE(1) + log(1 + exp(lngE(end) - lngE(1))) - log(4);
else
    lgC = lngE(end) + log(1 + exp(lngE(1) - lngE(end))) - log(4);
end
lngE = lngE - lgC;
for i = 1:length(lngE)
    if lngE(i) < 0
        lngE(i) = 0;
    end
end
        

figure(1)
plot(Energies, lngE)
figure(2)
plot(Energies, Hist)

function spins = randS(L)
spins = zeros(L);
for i = 1:L
    for j = 1:L
        spins(i, j) = sign(2 * rand() - 1);
    end
end
end

function E = CEnergy(spins, J, L)
% Energy of a 2D Ising lattice
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


