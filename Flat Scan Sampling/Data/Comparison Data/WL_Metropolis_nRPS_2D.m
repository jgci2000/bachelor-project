clear
close

% Comparison between WL and Metropolis and new RPS for 2D 4L lattice
% João Inácio, 18 Oct., 2020

load("JDOS_2D_4L_WL.mat", "g");
gCPP = load("JDOS_2D_4L_10E4_nRPS_CPP.txt");

L = 4;
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

gCPP(:, (NSpins / 2) + 2:end) = gCPP(:, (NSpins / 2):-1:1);

gECPP = zeros(1, NE);
for i = 1:NE
    gECPP(1, i) = sum(gCPP(i, :));
end

% Temperatures
T = 0.75:0.1:8;
NT = length(T);
kB = 1;
beta = zeros(1, NT);

Z = zeros(1, NT);
M = zeros(1, NT);

ZCPP = zeros(1, NT);
MCPP = zeros(1, NT);

for i = 1:NT
    beta(i) = - 1 / (T(i) * kB);
    
    % Partition Function
    Z(i) = sum(gE .* exp(energies * beta(i)));
    ZCPP(i) = sum(gECPP .* exp(energies * beta(i)));
    
    
    for idxE = 1:NE
        for idxM = 1:NM
            M(i) = M(i) + (g(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / Z(i);
            MCPP(i) = MCPP(i) + (gCPP(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / ZCPP(i);
        end
    end
    
end

figure(1)
% plot(T, Mmean/NSpins_Met)
hold on
plot(T, M/NSpins)
plot(T, MCPP/NSpins)
legend("WL", "new RPS(10E4 REP)")
title("Magnetization vs Temperature for 4L")
xlabel("Temperature")
ylabel("Magnetization")



