clear
close
clc

        % Wang-Landau Sampling for the Ising Model (DOS)

        % Defining some constants
% We need to select a starting value of f, the modification factor that it
% is not too large or too small, becuase it will lead to statistical erros
% or a long computation time, respectively. So f0 = e is a safe bet.
f = exp(1);     % Modification factor
n = 2;          % Exponent for the reduction of the modification factor.
L = 8;          % Number of spins for each side of the lattice.
N = L^2;        % Number of spins.
nConfig = 2^N;  % Number of different spin configurations
J = 1;          % Interaction strenght between particles.
% We need to stop the algorithm when all of the states have been visited a
% similar number of times, in other words, when the histogram is flat. As
% it is impossible to obtain a perfectly flat histogram, we check its
% "flatness" when H(E) is not less than p * average(H(E)).
p = 0.80;

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
EMin = computeEnergy(SMin, L, J);
EMax = - EMin;

% The energy is always a multiple of J. So we can create this vector of
% energies.
E = EMin:4*J:EMax;
NE = length(E);

% The density of states (DOS, g(E)) function is unknown at the begining of
% the program. It is our goal to find. The simplest approach to the proble 
% is to consider g(E) = 1.
g = ones(1, NE);
% Because g(E) becomes very large, very soon, we must work with the
% ln(g(E)). So...
lng = log(g);

% If the modification factor is less than exp(10^(-8)), stop. Else
% continue.
tic
while f > 1+1E-8
    % For each new value of f, we must reset the histogram.
    H = zeros(1, NE);
    
    % At the end of 10000 MCSweeps we must check if the histogram is flat.
    ex = 0;         % We do not wish to exit now!
    aux = 0;        % Initialize the aux variable.
    
    % If aux is empty, then the histogram is flat.
    while ~isempty(aux)
        % Each time this "for" loop runs, it executes one MC Sweep.
        for MCSweep = 1:10000
            for ni = 1:N
                % Now we have to choose a random spin and flip it, to calculate
                % the new energy of the system.
                i = randi([1 L]);
                j = randi([1 L]);
                SFlipped = S;
                SFlipped(i, j) = - SFlipped(i, j);
                
                % Compute the energy of S and SFlipped.
                E1 = computeEnergy(S, L, J);
                E2 = computeEnergy(SFlipped, L, J);
                
                % Index of the two energies
                idxE1 = findIndex(E, E1);
                idxE2 = findIndex(E, E2);
                
                % Probability of change
                % If P == 1, then SFlipped becomes S, or if P != 1, there 
                % is a probability (g(idxE1)/g(idxE2)) for the system to 
                % change. Becuase we are working with the ln(g), this 
                % g(idxE1)/g(idxE2) becomes this exp(lng(idxE1)-lng(idxE2)).
                ratio = exp(lng(idxE1) - lng(idxE2));
                P = min(ratio, 1);
                
                % When changing the system configuration, we also need to
                % update g and the histogram.
                if P == 1
                    H(idxE2) = H(idxE2) + 1;
                    S = SFlipped;
                    lng(idxE2) = lng(idxE2) + log(f);
                elseif P < 1
                    % We select E2 if, a random number between 0 and 1 is 
                    % less or equal than exp(lng(idxE1)-lng(idxE2)), else
                    %  we select E1.
                    if rand(1) <= ratio
                        H(idxE2) = H(idxE2) + 1;
                        S = SFlipped;
                        lng(idxE2) = lng(idxE2) + log(f);
                    else
                        H(idxE1) = H(idxE1) + 1;
                        lng(idxE1) = lng(idxE1) + log(f);
                    end
                end
            end
        end
        % Visualize the spins
%         figure(1)
%         imagesc(S)
%         drawnow
        
        % Visualize the DOS
        figure(2)
        lngT = lng - log(exp(lng(E == -2 * N))) + log(nConfig);
        plot(E/N, lngT, 'b')
        hold on
        plot(E/N, lng, 'r')
        pause(0.3)
        hold off
        
        % Compute the average of H(E)
        avgH = mean(H(H > 0));
        % As stated before, the minimum value of H(E) must be
        minH = avgH * p;
        aux = H(H(H > 0) < minH);
        
        % Plot the histogram
%         figure(2)
%         bar(E, H)
%         xlabel('E')
%         ylabel('Number of times visited')
%         pause(0.1)
    end
    
    % Update f
    f = f^(1/n);
end
toc

g = exp(lng);
% After this we can compute the partition function Z and with that, all
% sorts of thermodynamic quantities.

figure(2)
plot(E/N, lng)














