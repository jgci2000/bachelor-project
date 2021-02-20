clear
close

% Comparison between WL and Metropolis sampling for 2D 8L lattice
% João Inácio, 12 Sep., 2020

load("JDOS_8L_WL.mat", "g");
load("meanM_8L_Metropolis", "Mmean");

L = 8;
NSpins = L^2;

% All of the possible energies
energies = 4 * (0:NSpins) - 2 * NSpins;
NE = length(energies);

% Magnetizations
magnetizations = - NSpins:2:NSpins;
NM = length(magnetizations);

% Get the DOS
gE = zeros(1, NE);
for i = 1:NE
    gE(1, i) = sum(g(i, :));
end

% Temperatures
T = 0.75:0.1:8;
NT = length(T);
kB = 1;
beta = zeros(1, NT);

Z = zeros(1, NT);
M = zeros(1, NT);

for i = 1:NT
    beta(i) = - 1 / (T(i) * kB);
    
    % Partition Function
    Z(i) = sum(gE .* exp(energies * beta(i)));
    
    
    for idxE = 1:NE
        for idxM = 1:NM
            M(i) = M(i) + (g(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / Z(i);
        end
    end
    
end

figure(1)
plot(T, Mmean/NSpins)
hold on
plot(T, M/NSpins)
legend("Metropolis", "WL")
title("Magnetization vs Temperature for 8L")
xlabel("Temperature")
ylabel("Magnetization")



