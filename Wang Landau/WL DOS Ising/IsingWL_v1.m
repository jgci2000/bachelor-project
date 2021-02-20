clear
close
%clc

        % Wang-Landau Sampling for the Ising Model (DOS)

        % Defining some constants
% We need to select a starting value of f, the modification factor that it
% is not too large or too small, becuase it will lead to statistical erros
% or a long computation time, respectively. So f0 = e is a safe bet.
f = exp(1);     % Modification factor
n = 2;          % Exponent for the reduction of the modification factor.
L = 2;          % Number of spins for each side of the lattice.
N = L * L;        % Number of spins.
nConfig = 2^N;  % Number of different spin configurations
J = 1;          % Interaction strenght between particles.
% We need to stop the algorithm when all of the states have been visited a
% similar number of times, in other words, when the histogram is flat. As
% it is impossible to obtain a perfectly flat histogram, we check its
% "flatness" when H(E) is not less than p * average(H(E)).
p = 0.95;

        % Generate Spins Matrix
% It can be random. We don´t change it during the algorithm.
op = 3;
S = randomSpins(L, op);

        % W-L Algorithm

% For each value of energy there is a number os time visited. To check
% if the energy as already been visited, we must sotre its values.
% We have to discretize the energy domain to get the size of the vectors H
% and g, and the index of each correspondent energy.
% For that, we have to get EMin (all up or all down) and EMax (-EMin). 
SMin = randomSpins(L, 1);
EMin = computeEnergy1(SMin, L, J);
EMax = -EMin;

% The energy is always a multiple of J. So we can create this vector of
% energies.
temp = EMin:4*J:EMax;
E(1) = temp(1);
E(2:length(temp) - 3) = temp(3:end-2);
E(length(temp) - 2) = temp(end);
NE = length(E);

% The density of states (DOS, g(E)) function is unknown at the begining of
% the program. It is our goal to find. The simplest approach to the proble 
% is to consider g(E) = 1.
g = ones(1, NE);
% Because g(E) becomes very large, very soon, we must work with the
% ln(g(E)). So...
lngE = log(g);

% If the modification factor is less than exp(10^(-8)), stop. Else
% continue.
tic
while f > 1+1E-8
    % For each new value of f, we must reset the histogram.
    H = zeros(1, NE);
    
    % At the end of 10000 MCSweeps we must check if the histogram is flat.
    MCSweeps = 0;   % Reset the sweeps counter.
    ex = 0;         % We do not wish to exit now!
    
    % If exit == 1 then the histogram is flat.
    while ex ~= 1
        % Each time this "for" loop runs, it executes one MC Sweep.
        
            % Now we have to choose a random spin and flip it, to calculate
            % the new energy of the system.
            i = randi([1 L]);
            j = randi([1 L]);
            SFlipped = S;
            SFlipped(i, j) = - SFlipped(i, j);
            
            % Compute the energy of S and SFlipped.
            E1 = computeEnergy1(S, L, J);
            E2 = computeEnergy1(SFlipped, L, J);
            
            % Index of the two energies
            idxE1 = findIndex(E, E1);
            idxE2 = findIndex(E, E2);
            
            % Probability of change
            % If P == 1, then SFlipped becomes S, or if P != 1, there is a
            % probability (g(idxE1)/g(idxE2)) for the system to change.
            % Becuase we are working with the ln(g), this g(idxE1)/g(idxE2)
            % becomes this exp(lng(idxE1)-lng(idxE2)).
            ratio = exp(lngE(idxE1)-lngE(idxE2));
            P = min(ratio, 1);
                        
            % When changing the system configuration, we also need to 
            % update g and the histogram.
            if P == 1
                H(idxE2) = H(idxE2) + 1;
                S = SFlipped;
                lngE(idxE2) = lngE(idxE2) + log(f);
            elseif P < 1
                % We select E2 if, a random number between 0 and 1 is less
                % or equal than P, else we select E1.
                randNum = rand(1);
                if randNum <= ratio
                    H(idxE2) = H(idxE2) + 1;
                    S = SFlipped;
                    lngE(idxE2) = lngE(idxE2) + log(f);
                else
                    H(idxE1) = H(idxE1) + 1;
                    lngE(idxE1) = lngE(idxE1) + log(f);
                end
            end
        
        % Visualize the spins
%         figure(1)
%         imagesc(S)
%         drawnow
        
        MCSweeps = MCSweeps + 1;
        
        if MCSweeps == 10000
            MCSweeps = 0;
            % Compute the average of H(E)
            avgH = mean(H(H > 0));
            % As stated before, the minimum value of H(E) must be
            minH = avgH * p;
            aux = H(H(H > 0) < minH);
            
            % Plot the histogram
%             figure(2)
%             bar(E, H)
%             xlabel('E')
%             ylabel('Number of times visited')
%             pause(0.1)
            
            % If the aux vector is empty the histogram is flat.
            if isempty(aux)
                ex = 1;
            end
        end
    end
    
    % Update f
    f = f^(1/n);
end
toc

g = exp(lngE);
% After this we can compute the partition function Z and with that, all
% sorts of thermodynamic quantities.

lngN = lngE - lngE(E == -2 * N) + log(2);

% Normalize the DOS
if lngE(end) < lngE(1)
    lgC = lngE(1) + log(1+ exp(lngE(end) - lngE(1))) - log(4);
else
    lgC = lngE(end) + log(1+ exp(lngE(1) - lngE(end))) - log(4);
end
lngE = lngE - lgC;
for i = 1:length(lngE)
    if lngE(i) < 0
        lngE(i) = 0;
    end
end

figure(2)
plot(E/N, lngE)
hold on
plot(E/N, lngN)













