clear
close

% Comparison between WL and Metropolis sampling for 3D 4L lattice
% João Inácio, 5 Oct., 2020

load("JDOS_3D_4L_WL.mat", "JDOS");
load("meanM_3D_4L_Metropolis.mat", "Mmean");
g_CPP = load("JDOS_3D_4L_WL_CPP.txt");

L = 4;
NN = 6;
N_SPINS = L ^ 3;
g = JDOS;

% All of the possible energies
energies = - (1 / 2) * NN * N_SPINS:4:(1 / 2) * NN * N_SPINS;
NE = length(energies);

% Magnetizations
magnetizations = - N_SPINS:2:N_SPINS;
NM = length(magnetizations);

% Get the DOS
gE = zeros(1, NE);
for i = 1:NE
    gE(i) = sum(g(i, :));
end

gE_CPP = zeros(1, NE);
for i = 1:NE
    gE_CPP(i) = sum(g_CPP(i, :));
end

% Temperatures
T = 1.5:0.25:7.5;
NT = length(T);
kB = 1;
beta = zeros(1, NT);

Z = zeros(1, NT);
M = zeros(1, NT);

Z_CPP = zeros(1, NT);
M_CPP = zeros(1, NT);

for i = 1:NT
    beta(i) = - 1 / (T(i) * kB);
    
    % Partition Function
    Z(i) = sum(gE .* exp(energies * beta(i)));
    Z_CPP(i) = sum(gE_CPP .* exp(energies * beta(i)));
    
    for idxE = 1:NE
        for idxM = 1:NM
            M(i) = M(i) + (g(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / Z(i);
            M_CPP(i) = M_CPP(i) + (g_CPP(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / Z_CPP(i);
        end
    end
    
end

figure(1)
plot(T, Mmean/N_SPINS)
hold on
plot(T, M/N_SPINS)
plot(T, M_CPP/N_SPINS)
legend("Metropolis", "WL MATLAB", "WL CPP")
title("Magnetization vs Temperature for 4L")
xlabel("Temperature")
ylabel("Magnetization")



