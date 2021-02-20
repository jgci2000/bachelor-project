clear
close
clc

    % Metropolis Algorithm for the Ising Model

    % Some variables
L = 32;          % Number of spins for each side of the lattice.
N = L^2;        % Number of spins.
nConfig = 2^N;  % Number of different spin configurations
J = 1;          % Interaction strenght between particles.
T = 2.6;          % Temperature. We can change this value to observe some
% magnetic properties.
kB = 1;         % Boltzman constant in reduce units.
beta = 1 / (T * kB);
itMax = 1E6;  % Maximum number of iterations for the method. 
% 100 MCSweeps.

% Generate the first spin configuration
% We get a L*L matrix with random +/-1 spins.
S = ones(L);

    % Metropolis Method
% Mean magnetization vector.
M = zeros(1, itMax);
% Energy per iteration vector.
H = zeros(1, itMax);

% We repeat the method a suficient number of times. This value may need to
% be ajusted.
tic
for it = 1:itMax
    % Flip a random spin.
    i = randi([1 L]);
    j = randi([1 L]);
    SFlipped = S;
    SFlipped(i, j) = -SFlipped(i, j);
    
    % Compute the new energy.
    E1 = CEnergy(S, J, L);
    E2 = CEnergy(SFlipped, J, L);
    
    % Difference of energies.
    dE = E2 - E1;
    
    % Random number r for the acception or rejection of the new
    % configuration.
    r = rand(1);
    ratio = exp(-dE * beta);
    
    if r < ratio
        S = SFlipped;
    end
    
    % Compute the mean magnetization for the kth iteration.
    M(it) = sum(sum(S));
end
toc

% Some Plots
figure(2)
plot(1:itMax, M/N)
xlabel('Number of Iterations')
ylabel('Magnetization')













