function g = WL_JDOS(L, fMin)

NSpins = L * L;     % Number of spins
J = 1;              % Interaction strenght between particles

    % Neighbours Vectors
aux = 1:L;
menos = circshift(aux, [0 -1]);
mais = circshift(aux, [0 +1]);

    % Wang Landau Sampling
f = exp(1);     % Modification factor
p = 0.90;       % Flatness -> min(hist) > avg(hist)*p

% All of the possible energies
energies = 4 * (0:NSpins) - 2 * NSpins;
NE = NSpins + 1;

% All of the possible magnetizations. 
magnetizations = - NSpins:2:NSpins;
NM = NSpins + 1;

% Random configuration
spins = randSpins(L);
E = CEnergy(spins, J, L);
M = CMagnetization(spins, L);
idxE = find(energies == E);
idxM = find(magnetizations == M);

% Joint Desnsity of States Vector
g = ones(NE, NM);
% Rows -> Energies columns -> Magnetizations
% Working with the ln(g(E)) is better because the values are too high.
lngEM = log(g);

% At the end of 10000 MCSweeps check if the histogram is flat.
MCSweeps = 0;   % Reset the sweeps counter.

% Histogram
hist = zeros(NE, NM);
% Rows -> Energies columns -> Magnetizations

% Start time
tic;

while f > fMin
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
        MNew = M - 2 * S;
        idxEN = find(energies == ENew);
        idxMN = find(magnetizations == MNew);
        
        % Flip the spin
        ratio = exp(lngEM(idxE, idxM) - lngEM(idxEN, idxMN));
        P = min(ratio, 1);
        
        if P > rand()
            spins(i, j) = -S;
            E = ENew;
            M = MNew;
            idxE = idxEN;
            idxM = idxMN;
        end
        
        % Update the histogram and g(E)
        hist(idxE, idxM) = hist(idxE, idxM) + 1;
        lngEM(idxE, idxM) = lngEM(idxE, idxM) + log(f);
    end
    
    MCSweeps = MCSweeps + 1;
    
    % Check the flatness of the histogram each 10000 setps.
    if mod(MCSweeps, 10000) == 0
        avgH = mean(mean(hist(hist > 0)));
        minH = min(min(hist(hist > 0)));
        
        if minH > avgH * p
            hist = zeros(NE, NM);
%             fprintf("%d: the histogram is flat. Min: %.0f Avg: %.2f f: %.8f\n", MCSweeps, minH, avgH, f);
            f = f^(1/2);
            % Reset the counter.
            MCSweeps = 0;
        end
    end
end

% End time
endT = toc;

% Normalize the DOS
% lngE is the normalized Density of States.
lngEM(lngEM > 0) = lngEM(lngEM > 0) - lngEM(energies == -2 * NSpins, magnetizations == - NSpins) + log(2);

% Get the actual JDOS
g = zeros(NE, NM);

for i = 1:NE
    for j = 1:NM
        if lngEM(i, j) ~= 0
            g(i, j) = exp(lngEM(i, j)) / 2;
        end
    end
end

% The runtime is the last entry in the JDOS matrix
g(1, NM + 1) = endT;

end

